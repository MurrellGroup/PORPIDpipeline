# load required packages
using NextGenSeqUtils, StatsBase
using BioSequences, FASTX
using CodecZlib: GzipDecompressorStream
using CodecZlib: GzipCompressorStream

function pp_demux_dict(names, seqs, fwd_primers, rev_primers, reject_trim_writer;
    verbose = true, phreds = nothing, tol_one_error = true, demux_dir = "demux")
    if rev_primers == nothing
        fwd_matches = fast_primer_match(seqs,fwd_primers,tol_one_error=tol_one_error)
        rev_comp_bool = fwd_matches .< 0
        keepers = abs.(fwd_matches) .> 0
        fwd_matches = abs.(fwd_matches)
        pair_keeps = fwd_matches[keepers]
    else
         keepers,fwd_matches,rev_matches,rev_comp_bool = fast_primer_pair_match(seqs,fwd_primers,rev_primers,tol_one_error=tol_one_error)
        f_keeps = fwd_matches[keepers]
        r_keeps = rev_matches[keepers]
        pair_keeps = [(f_keeps[i],r_keeps[i]) for i in 1:length(f_keeps)]
    end
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            println(s[1], " => ", s[2])
        end
    end

    if phreds == nothing
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                if rev_primers == nothing
                    d_key = fwd_matches[i]
                else
                    d_key = (fwd_matches[i],rev_matches[i])
                end
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            else
                # write_fasta(demux_dir*"/REJECTS_PRIMER_TRIM.fasta",
                #   [seqs[i]],names=["reject-$(i)"],append=true)
                # record = FASTA.Record("reject-$(i)", seqs[i])
                record = FASTA.Record("$(names[i])-rejected", seqs[i])
                write(reject_trim_writer, record)
            end
        end
        return seq_dict
    else
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Vector{Int8},Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                if rev_primers == nothing
                    d_key = fwd_matches[i]
                else
                    d_key = (fwd_matches[i],rev_matches[i])
                end
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            else
                # write_fastq(demux_dir*"/REJECTS_PRIMER_TRIM.fastq",
                #   [seqs[i]],[phreds[i]],names=["reject-$(i)"],append=true)
                # record = FASTQ.Record("reject-$(i)", seqs[i], phreds[i])
                record = FASTQ.Record("$(names[i])-rejected", seqs[i], phreds[i])
                write(reject_trim_writer, record)
            end
        end
        return seq_dict
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

function iterative_primer_match(seqs,full_primers,
  window::Int,slide_by::Int;tol_one_error=true)
    if(slide_by + window - 1 > minimum(length.(full_primers)))
        @warn("Matching window extends beyond shortest primer. This is ok, but check that you aren't matching something too short.",maxlog=1)
    end
    primers = [p[1:min(window,minimum(length.(full_primers)))] for p in full_primers]
    filter = fast_primer_match(seqs,primers,tol_one_error=tol_one_error);
    for i in 2:slide_by
        unresolved = filter .== 0
        primers = [p[i:min(i + window - 1,minimum(length.(full_primers)))] for p in full_primers]
        filter[unresolved] = fast_primer_match(seqs[unresolved],primers,tol_one_error=tol_one_error);
    end
    return filter
end


function sliding_demux_dict(names,seqs,fwd_primers,window::Int,slide_by::Int,reject_fwd_writer;
        verbose = true, phreds = nothing, tol_one_error = true, demux_dir = "demux")
    fwd_matches = iterative_primer_match(seqs,fwd_primers,window,slide_by,tol_one_error=tol_one_error)
    rev_comp_bool = fwd_matches .< 0
    keepers = abs.(fwd_matches) .> 0
    fwd_matches = abs.(fwd_matches)
    pair_keeps = fwd_matches[keepers]
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            println(s[1], " => ", s[2])
        end
    end
    if phreds == nothing
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]

                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],i))
                end
            else
                # write_fasta(demux_dir*"/REJECTS_PRIMER_FWD.fasta",
                #   [seqs[i]],names=["reject-$(i)"],append=true)
                # record = FASTA.Record("reject-$(i)", seqs[i])
                record = FASTA.Record("$(names[i])-rejected", seqs[i])
                write(reject_fwd_writer, record)
            end
        end
        return seq_dict
    else
        seq_dict = Dict()
        for pair in sorted_pairs
            seq_dict[pair[1]] = Tuple{String,Vector{Int8},Int64}[]
        end
        for i in 1:length(keepers)
            if keepers[i]
                d_key = fwd_matches[i]
                if rev_comp_bool[i]
                    push!(seq_dict[d_key],(reverse_complement(seqs[i]),reverse(phreds[i]),i))
                else
                    push!(seq_dict[d_key],(seqs[i],phreds[i],i))
                end
            else
                # write_fastq(demux_dir*"/REJECTS_PRIMER_FWD.fastq",
                #   [seqs[i]],[phreds[i]],names=["reject-$(i)"],append=true)
                # record = FASTQ.Record("reject-$(i)", seqs[i], phreds[i])
                record = FASTQ.Record("$(names[i])-rejected", seqs[i], phreds[i])
                write(reject_fwd_writer, record)
            end
        end
        return seq_dict
    end
end


function longest_conserved_5p(seqs)
    for i in 1:length(seqs[1])
        if length(unique(getindex.(seqs,i))) != 1
            return seqs[1][1:i-1]
        end
    end
    return seqs[1]
end

function chunked_filter_apply(in_path, out_path, func::Function; chunk_size=10000, f_kwargs = [], verbose = false)
    if endswith(in_path, ".gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(in_path)))
        reject_qual_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_DEMUX_QUAL.fastq.gz","w")))
        reject_rev_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_PRIMER_REV.fastq.gz","w")))
        reject_fwd_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_PRIMER_FWD.fastq.gz","w")))
        reject_trim_writer = FASTQ.Writer(GzipCompressorStream(open(out_path*"/REJECTS_PRIMER_TRIM.fastq.gz","w")))
    else
        reader = FASTQ.Reader(open(in_path))
        reject_qual_writer = FASTQ.Writer(open(out_path*"/REJECTS_DEMUX_QUAL.fastq","w"))
        reject_rev_writer = FASTQ.Writer(open(out_path*"/REJECTS_PRIMER_REV.fastq","w"))
        reject_fwd_writer = FASTQ.Writer(open(out_path*"/REJECTS_PRIMER_FWD.fastq","w"))
        reject_trim_writer = FASTQ.Writer(open(out_path*"/REJECTS_PRIMER_TRIM.fastq","w"))
    end
    # if !hasmethod(func, Int64, Int64, Tuple{Array{Any,1}, Array{Phred,1}, Array{Any,1}})
    #    @error "Function argument must accept func(seqs::Array{Any,1}, phreds::Array{Phred,1}, names::Array{Any,1})!"
    # end
    seqs, phreds, names = [], Vector{Phred}[], []
    i = 0
    read_counts = [0, 0, 0, 0, 0, 0]
    chunk = 0
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        push!(seqs, FASTQ.sequence(String, record))
        push!(phreds, FASTQ.quality(record, :sanger))
        push!(names, FASTQ.identifier(record))
        i += 1
        if i == chunk_size
            #apply func...
            chunk += 1
            println("processing chunk $(chunk), of size $(chunk_size)")
            read_counts += func(chunk, chunk_size, seqs, phreds, names,
                reject_qual_writer, reject_rev_writer, reject_fwd_writer, reject_trim_writer;
                f_kwargs...)
            seqs, phreds, names = [], Vector{Phred}[], []
            i = 0
        end
    end
    if i > 0
        #apply func...
        chunk += 1
        println("processing chunk $(chunk), of size $(i)")
        read_counts += func(chunk, chunk_size, seqs, phreds, names,
            reject_qual_writer, reject_rev_writer, reject_fwd_writer, reject_trim_writer;
            f_kwargs...)
    end
    close(reader)
    close(reject_qual_writer)
    close(reject_rev_writer)
    close(reject_fwd_writer)
    close(reject_trim_writer)
    return read_counts   # [total_reads, quality_reads, bad_reads, short_reads, long_reads, demuxed_reads]
end

#------Chunked quality demux function--------
"""
function chunked_quality_demux(chunk, chunk_size, seqs, phreds, names,
    reject_qual_writer, reject_rev_writer, reject_fwd_writer, reject_trim_writer;
    demux_dir = "demux", samples = Dic(), error_rate = 0.01, min_length = 30,
    max_length = 1000000, label_prefix = "seq", error_out = true, verbose = false)

Intended to be passed to `chunked_filter_apply()` for quality filtering and demultiplexing of
FASTQ files in one pass.
"""
function chunked_quality_demux(chunk, chunk_size, seqs, phreds, names,
    reject_qual_writer, reject_rev_writer, reject_fwd_writer, reject_trim_writer;
    demux_dir = "demux",
    samples = Dic(),
    error_rate = 0.01,
    min_length = 30,
    max_length = 1000000,
    label_prefix = "seq",
    error_out = true,
    verbose = true)
    
    # println("julia 1.7")

    total_reads, quality_reads, bad_reads, short_reads, long_reads, demuxed_reads = 0, 0, 0, 0, 0, 0
    
    total_reads = length(seqs)

    # quality filter
    if verbose
        println("filtering chunk on mean phred scores ...")
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
    
    # rename records using increasing seq no and annotaion
    if error_out == true
        names = ["$label_prefix$((chunk-1)*chunk_size+i)|ee=$(mean_errors[i])" for i in 1:length(names)]
    else
        names = ["$label_prefix$((chunk-1)*chunk_size+i)|" for i in 1:length(names)]
    end
    
    # save the rejected sequences as a fastq file
    rejects = setdiff(1:length(seqs),inds)
    # write_fastq(demux_dir*"/REJECTS_DEMUX_QUAL.fastq",
    #     seqs[rejects],phreds[rejects],names=names[rejects],append=true)
    for i in rejects
        record = FASTQ.Record(names[i], seqs[i], phreds[i])
        write(reject_qual_writer, record)
    end
    # process the passed sequences
    seqs, phreds, names = seqs[inds], phreds[inds], names[inds]
    
    quality_reads = length(seqs)
    
    #demux...
    if verbose
        println("de-multiplexing chunk ...")
    end
    #forward
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
        println("Splitting by primers...")
    end
    rev_adapters = [v["cDNA_primer"] for (k,v) in samples]
    rev_adapter = longest_conserved_5p(rev_adapters)
    rev_adapter = String(split(rev_adapter,r"[a-z]")[1]) #if all contain sample ID keep sample ID
    #sliding window demultiplex on forward primers
    fwd_demux_dic = sliding_demux_dict(names, seqs,
                                       unique_fwd_ends,
                                       12,
                                       13, #must not extend into the sampleID
                                       reject_fwd_writer,
                                       verbose=false,
                                       phreds = phreds,
                                       demux_dir = demux_dir)
    #iterate by each forward primer group
    for j in unique(fwd_end_group_arr)
        #define templates
        template_names = [k for (k,v) in samples][fwd_end_group_arr .== j]
        templates = [samples[n]["cDNA_primer"] for n in template_names]
        sampleIDs = uppercase.([m.match for m in match.(r"[a-z]+", templates)])
        IDind2name = Dict(zip(collect(1:length(sampleIDs)),template_names));

        #retrieve from demux_dic
        
        # hack to prevent missing key error
        if j in keys(fwd_demux_dic)   # start of hack
        
        seqs_fwd = [i[1] for i in fwd_demux_dic[j]];
        phreds_fwd = [i[2] for i in fwd_demux_dic[j]];
        seq_names_fwd = names[[i[3] for i in fwd_demux_dic[j]]]

        if verbose
            println("$(length(seqs_fwd)) reads matching forward primer $(unique_fwd_ends[j])")
        end

        #match to reverse adapter
        rev_matches = iterative_primer_match(seqs_fwd, [rev_adapter], 12, 13, tol_one_error=true);
        rev_keepers = rev_matches .< 0
        
        # save the reverse primer rejects as a fastq file
        rejects = .!(rev_keepers)
        # write_fastq(demux_dir*"/REJECTS_PRIMER_REV.fastq",
        #     seqs_fwd[rejects],phreds_fwd[rejects],names=seq_names_fwd[rejects],append=true)
        for i in 1:length(rejects)
            if rejects[i]
                record = FASTQ.Record("$(seq_names_fwd[i])-rejected", seqs_fwd[i], phreds_fwd[i])
                write(reject_rev_writer, record)
            end
        end
        #filter to reverse adapter matches
        seqs_both = seqs_fwd[rev_keepers]
        phreds_both = phreds_fwd[rev_keepers]
        seq_names_both = seq_names_fwd[rev_keepers]
        if verbose
            println("$(length(seqs_both)) contain reverse adapter $(rev_adapter)")
        end

        #trim and orient primers, phreds
        if verbose
            println("Trimming primers for read group $(j)...")
        end
        trimmed = [double_primer_trim(seqs_both[i],
                                            phreds_both[i],
                                            uppercase(unique_fwd_ends[j]),
                                            uppercase(rev_adapter)) for i in 1:length(seqs_both)]

        #separate by donor (sample) IDs
        if verbose
            println("Splitting read group $(j) by donor ID...")
        end
        split_donor_dic = pp_demux_dict(seq_names_both,
                                           [s[1] for s in trimmed],
                                           sampleIDs,
                                           nothing,
                                           reject_trim_writer,
                                           verbose=false,
                                           phreds=[s[2] for s in trimmed],
                                           tol_one_error=false,
                                           demux_dir = demux_dir)

        if verbose
            println("Writing individual donor files...")
        end
        for i in 1:length(sampleIDs)
            if i in keys(split_donor_dic)
                if verbose
                    println(IDind2name[i]," => ",length(split_donor_dic[i]))
                end
                template = template_names[i]
                write_fastq(demux_dir*"/"*template*".fastq",
                            [i[1] for i in split_donor_dic[i]],
                            [i[2] for i in split_donor_dic[i]];
                            names = seq_names_both[[i[3] for i in split_donor_dic[i]]],
                            append=true)
                demuxed_reads += length(split_donor_dic[i])
            else
                if verbose
                    println(IDind2name[i]," => 0")
                end
            end
        end
        end   # end of hack
    end
    return [total_reads, quality_reads, bad_reads, short_reads, long_reads, demuxed_reads]

end
