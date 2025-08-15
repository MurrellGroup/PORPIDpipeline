# load required packages
using BioSequences, FASTX, StatsBase
using CodecZlib: GzipDecompressorStream
using CodecZlib: GzipCompressorStream

function findbestlast(query,tol,window)
    best=nothing
    next=findfirst(query,tol,window)
    while !isnothing(next)
        best=next
        next=findnext(query,tol,window,last(next)+1)
    end
    if isnothing(best)
        return 0
    else
        return last(best) + 1
    end
end


function unique_not_substr(a)
    out = []
    for i in unique(a)
        res = true
        for j in unique(a)
            if occursin(i, j) & (i != j)
                res = false
            end
        end
        if res
            push!(out, i)
        end
    end
    return out
end

function longest_conserved_5p(seqs)
    for i in 1:length(seqs[1])
        if length(unique(getindex.(seqs,i))) != 1
            return seqs[1][1:i-1]
        end
    end
    return seqs[1]
end

function chunked_filter_apply(in_path, out_path, func::Function; chunk_size=10000, f_kwargs = [])
    if endswith(in_path, ".gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(in_path)))
        reject_qual_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_QUAL.fastq.gz","w")))
        reject_primer_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_PRIMER.fastq.gz","w")))
        reject_idtrim_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_IDTRIM.fastq.gz","w")))
    else
        reader = FASTQ.Reader(open(in_path))
        reject_qual_writer = FASTQ.Writer(open(out_path*"/REJECTS_QUAL.fastq","w"))
        reject_primer_writer = FASTQ.Writer(open(out_path*"/REJECTS_PRIMER.fastq","w"))
        reject_idtrim_writer = FASTQ.Writer(open(out_path*"/REJECTS_IDTRIM.fastq","w"))
    end
    
    seqs, phreds, names = [], Vector{Phred}[], []
    i = 0
    read_counts = [0, 0, 0, 0, 0, 0]
    chunk = 0
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        push!(seqs, FASTQ.sequence(LongDNA{4}, record))
        push!(phreds, collect(FASTQ.quality_scores(record, :sanger)))  # quality changed to quality_scores
        push!(names, FASTQ.identifier(record))
        i += 1
        if i == chunk_size
            #apply func...
            chunk += 1
            println("processing chunk $(chunk), of size $(chunk_size)")
            read_counts += func(chunk, chunk_size, seqs, phreds, names,
                reject_qual_writer, reject_primer_writer, reject_idtrim_writer;
                f_kwargs...)
            @show read_counts
            seqs, phreds, names = [], Vector{Phred}[], []
            i = 0
        end
    end
    if i > 0
        #apply func...
        chunk += 1
        println("processing chunk $(chunk), of size $(i)")
        read_counts += func(chunk, chunk_size, seqs, phreds, names,
            reject_qual_writer, reject_primer_writer, reject_idtrim_writer;
            f_kwargs...)
        @show read_counts
    end
    close(reader)
    close(reject_qual_writer)
    close(reject_primer_writer)
    close(reject_idtrim_writer)
    return read_counts   # [total_reads, quality_reads, bad_reads, short_reads, long_reads, demuxed_reads]
end

#------Chunked quality demux function--------

function chunked_quality_demux(chunk, chunk_size, seqs, phreds, names,
    reject_qual_writer, reject_primer_writer, reject_idtrim_writer;
    demux_dir = "demux",
    samples = Dic(),
    error_rate = 0.01,
    min_length = 30,
    max_length = 1000000,
    primer_tol=1,
    primer_window=100,
    primer_chop=0,
    label_prefix = "seq",
    error_out = true,
    verbose = false)

    total_reads, quality_reads, bad_reads, short_reads, long_reads, demuxed_reads = 0, 0, 0, 0, 0, 0
    
    total_reads = length(seqs)

    # quality filter
    if verbose
        println("filtering chunk of size $(total_reads) on mean phred scores ...")
    end
    lengths = length.(seqs)
    mean_errors = [mean(phred_to_p.(phred)) for phred in phreds]
    bad_inds = mean_errors .>= error_rate
    good_inds = mean_errors .< error_rate
    short_inds = lengths .<= min_length
    long_inds = lengths .>= max_length
    bad_reads = sum(bad_inds)
    short_reads = sum( (short_inds) .& (good_inds) )
    long_reads = sum( (long_inds) .& (good_inds) )
    inds = [1:length(seqs);][(lengths .< max_length) .& (lengths .> min_length) .& (mean_errors .< error_rate)]
    good_reads = length( inds )
    
    if verbose
        @show good_reads, bad_reads, short_reads, long_reads
    end
    
    # rename records using increasing seq no and annotaion
    if error_out == true
        names = ["$label_prefix$((chunk-1)*chunk_size+i)|ee=$(mean_errors[i])" for i in 1:length(names)]
    else
        names = ["$label_prefix$((chunk-1)*chunk_size+i)|" for i in 1:length(names)]
    end
    
    # save the rejected sequences as a fastq file
    rejects = setdiff(1:length(seqs),inds)
    for i in rejects
        record = FASTQ.Record(names[i], seqs[i], phreds[i])
        write(reject_qual_writer, record)
    end
    # process the passed sequences
    seqs, phreds, names = seqs[inds], phreds[inds], names[inds]
    
    quality_reads = length(seqs)
    
    #demux...
    if verbose
        println("de-multiplexing chunk with $(length(seqs)) quality reads ...")
    end
    
    fwd_ends = [v["sec_str_primer"] for (k,v) in samples]
    unique_fwd_ends = unique_not_substr(fwd_ends)
    fwd_end_group_arr = []
    for e in fwd_ends
        matches = []
        for (j, group) in enumerate(unique_fwd_ends)
            if occursin(e, group)
                push!(matches, (group => j))
            end
        end
        push!(fwd_end_group_arr, sort(matches, by = x -> length(x[1]))[1][2])
    end
    if verbose
        @show fwd_ends, unique_fwd_ends, fwd_end_group_arr
    end

    rev_adapters = [v["cDNA_primer"] for (k,v) in samples]
    rev_adapter = longest_conserved_5p(rev_adapters)
    rev_adapter = String(split(rev_adapter,r"[a-z]")[1]) #if all contain sample ID keep sample ID
    if verbose
        println("Splitting by primers...")
        @show unique_fwd_ends
        @show rev_adapter
    end
    
    ###### new code using ApproximateSearchQuery from BioSequences
    window=primer_window
    tol=primer_tol
    chop=primer_chop
    
    if verbose
        println("searching for primers chopped by $(chop) nucs, using tol=$(tol) and window=$(window)")
    end
    
    fwd_rejects=0
    rev_rejects=0
    demuxed_reads=0
    all_keeps=[]
    for j in unique(fwd_end_group_arr)
        #define templates
        template_names = [k for (k,v) in samples][fwd_end_group_arr .== j]
        templates = [samples[n]["cDNA_primer"] for n in template_names]
        sampleIDs = uppercase.([m.match for m in match.(r"[a-z]+", templates)])
        fwd_primers = [samples[n]["sec_str_primer"] for n in template_names]
        fwd_primer = union(fwd_primers)[1]
    
        IDind2name = Dict(zip(collect(1:length(sampleIDs)),template_names));
        if verbose
            @show sampleIDs, IDind2name, fwd_primers, fwd_primer
        end
        
        # fwd_queries=[ApproximateSearchQuery(LongDNA{4}(fwd)) for fwd in fwd_ends]
        fwd_query=ApproximateSearchQuery(LongDNA{4}(fwd_primer[chop+1:end]))
        rev_query=ApproximateSearchQuery(LongDNA{4}(rev_adapter[chop+1:end]))
            
        keeps=[]
        for i in 1:length(seqs)
            rev=false
            fwd_window=seqs[i][1:window]
            fwd_trim = findbestlast(fwd_query, tol, fwd_window)
            # no luck, try reverse complement
            if fwd_trim==0
                BioSequences.reverse_complement!(seqs[i])
                reverse!(phreds[i])
                fwd_window=seqs[i][1:window]
                # fwd_trim = maximum((x->findbestlast(x, tol, fwd_window )).(fwd_queries))
                fwd_trim = findbestlast(fwd_query, tol, fwd_window)
                rev=true
            end
            if fwd_trim==0
                fwd_rejects+=1
                # record = FASTQ.Record("$(names[i])-nofwdprimer", seqs[i], phreds[i])
                # write(reject_fwd_writer, record)
            else # look for the reverse primer
                rev_window=seqs[i][end-window+1:end]
                rev_trim=findbestlast(rev_query, tol, BioSequences.reverse_complement(rev_window))
                if rev_trim==0
                    rev_rejects+=1
                    # record = FASTQ.Record("$(names[i])-norevprimer", seqs[i], phreds[i])
                    # write(reject_rev_writer, record)
                else
                    # primers present, now trim primers from sequence and phreds
                    seqs[i]=seqs[i][fwd_trim:end-rev_trim+1]
                    phreds[i]=phreds[i][fwd_trim:end-rev_trim+1]
                    push!(keeps,i)
                    push!(all_keeps,i)
                end
            end
        end
        seqs_keep, phreds_keep, names_keep = BioSequences.reverse_complement.(seqs[keeps]), reverse.(phreds[keeps]), names[keeps]
        #separate by donor (sample) IDs
        if verbose
            println("Splitting read group $(j) by donor ID...")
            println("Writing individual donor files...")
        end
        
        id_keeps=[]
        for i in 1:length(sampleIDs)
            template = template_names[i]
            stream = FASTQ.Writer(GzipCompressorStream(open(demux_dir*"/"*template*".fastq.gz", append=true)))
            id_query=ApproximateSearchQuery(LongDNA{4}(sampleIDs[i]))
            count=0
            for k in 1:length(seqs_keep)
                id_hit=findfirst(id_query,0,seqs_keep[k][1:10])
                if !isnothing(id_hit)
                    start=first(id_hit)
                    write(stream, FASTQ.Record(names_keep[k], seqs_keep[k][start:end], phreds_keep[k][start:end]))
                    count+=1
                    demuxed_reads+=1
                    push!(id_keeps,k)
                end
            end
            @show template, count
            close(stream)
        end
        # save idtrim rejects
        id_rejects = setdiff(1:length(seqs_keep),id_keeps)
        for k in id_rejects
            record = FASTQ.Record("$(names_keep[k])-noid", seqs_keep[k], phreds_keep[k])
            write(reject_idtrim_writer, record)
        end
        if verbose
            println("primer group $(j)...")
            @show length(id_rejects)
        end
    end
    
    # save the rejected sequences
    not_keeps= collect( setdiff(Set(1:length(seqs)),Set(all_keeps)) )
    for i in not_keeps
        record = FASTQ.Record("$(names[i])-noprimerpair", seqs[i], phreds[i])
        write(reject_primer_writer, record)
    end
    
    if verbose
        @show length(not_keeps), demuxed_reads
    end
    
    return [total_reads, quality_reads, bad_reads, short_reads, long_reads, demuxed_reads]
    
end

