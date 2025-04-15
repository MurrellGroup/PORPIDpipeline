# load required packages
using BioSequences, FASTX, StatsBase
using CodecZlib: GzipDecompressorStream
using CodecZlib: GzipCompressorStream

function pp_demux_dict(names, seqs, fwd_primers, rev_primers, reject_trim_writer;
    verbose = false, phreds = nothing, tol_one_error = true, demux_dir = "demux")
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
            @show s
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
  window::Int,slide_by::Int;tol_one_error=true,verbose=false)
    if(slide_by + window - 1 > minimum(length.(full_primers)))
        @warn("Matching window extends beyond shortest primer. This is ok, but check that you aren't matching something too short.",maxlog=1)
    end
    primers = [p[1:min(window,minimum(length.(full_primers)))] for p in full_primers]
    if verbose
        @show full_primers, primers
    end
    filter = fast_primer_match(seqs,primers,tol_one_error=tol_one_error);
    for i in 2:slide_by
        unresolved = filter .== 0
        primers = [p[i:min(i + window - 1,minimum(length.(full_primers)))] for p in full_primers]
        filter[unresolved] = fast_primer_match(seqs[unresolved],primers,tol_one_error=tol_one_error);
    end
    return filter
end


function sliding_demux_dict(names,seqs,fwd_primers,window::Int,slide_by::Int,reject_fwd_writer;
        verbose = false, phreds = nothing, tol_one_error = true, demux_dir = "demux")
    fwd_matches = iterative_primer_match(seqs,fwd_primers,window,slide_by,tol_one_error=tol_one_error)
    rev_comp_bool = fwd_matches .< 0
    keepers = abs.(fwd_matches) .> 0
    fwd_matches = abs.(fwd_matches)
    pair_keeps = fwd_matches[keepers]
    pair_counts = countmap(pair_keeps)
    sorted_pairs = sort([(k,pair_counts[k]) for k in keys(pair_counts)])
    if verbose
        for s in sorted_pairs
            @show s
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

function chunked_filter_apply(in_path, out_path, func::Function; chunk_size=10000, f_kwargs = [])
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
        push!(phreds, collect(FASTQ.quality_scores(record, :sanger)))  # quality changed to quality_scores
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
    verbose = false)
    
    # println("julia 1.7")

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
    if verbose
        @show rev_adapter
    end
    #sliding window demultiplex on forward primers
    fwd_demux_dic = sliding_demux_dict(names, seqs,
                                       unique_fwd_ends,
                                       12,
                                       13, #must not extend into the sampleID
                                       reject_fwd_writer,
                                       verbose=verbose,
                                       phreds = phreds,
                                       demux_dir = demux_dir)
    #iterate by each forward primer group
    if verbose
        @show keys(fwd_demux_dic)
        @show unique_fwd_ends
        @show unique(fwd_end_group_arr)
    end
    for j in unique(fwd_end_group_arr)
        #define templates
        template_names = [k for (k,v) in samples][fwd_end_group_arr .== j]
        templates = [samples[n]["cDNA_primer"] for n in template_names]
        sampleIDs = uppercase.([m.match for m in match.(r"[a-z]+", templates)])
        if verbose
            @show sampleIDs
        end

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
        rev_matches = iterative_primer_match(seqs_fwd, [rev_adapter], 12, 13, tol_one_error=true, verbose=verbose);
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
                                           verbose=verbose,
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
 #                write_fastq(demux_dir*"/"*template*".fastq",
 #                            [i[1] for i in split_donor_dic[i]],
 #                            [i[2] for i in split_donor_dic[i]];
 #                            names = seq_names_both[[i[3] for i in split_donor_dic[i]]],
 #                            append=true)
                stream = FASTQ.Writer(GzipCompressorStream(open(demux_dir*"/"*template*".fastq.gz", append=true)))
                for (s, q, n) in zip([i[1] for i in split_donor_dic[i]],
                                     [i[2] for i in split_donor_dic[i]],
                                     seq_names_both[[i[3] for i in split_donor_dic[i]]])
                    write(stream, FASTQ.Record(n, s, q))
                end
                close(stream)
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

# some more functions copypasta from demux.jl of NextGenSeqUtils

const IUPACnumDict = Dict{Char,Int8}(Dict())
IUPACnumDict['A']=1
IUPACnumDict['C']=2
IUPACnumDict['G']=3
IUPACnumDict['T']=4
IUPACnumDict['U']=5
IUPACnumDict['R']=6
IUPACnumDict['Y']=7
IUPACnumDict['S']=8
IUPACnumDict['W']=9
IUPACnumDict['K']=10
IUPACnumDict['M']=11
IUPACnumDict['B']=12
IUPACnumDict['D']=13
IUPACnumDict['H']=14
IUPACnumDict['V']=15
IUPACnumDict['N']=16

IUPACnucLookup=[[true false false false false true false false true false true false true true true true]
[false true false false false false true true false false true true false true true true]
[false false true false false true false true false true false true true false true true]
[false false false true true false true false true true false true true true false true]
[false false false true true false true false true true false true true true false true]
[true false true false false true false true true true true true true true true true]
[false true false true true false true true true true true true true true true true]
[false true true false false true true true false true true true true true true true]
[true false false true true true true false true true true true true true true true]
[false false true true true true true true true true false true true true true true]
[true true false false false true true true true false true true true true true true]
[false true true true true true true true true true true true true true true true]
[true false true true true true true true true true true true true true true true]
[true true false true true true true true true true true true true true true true]
[true true true false false true true true true true true true true true true true]
[true true true true true true true true true true true true true true true true]];

function toIUPACnum(strin::String)
    return [IUPACnumDict[i] for i in strin]
end

function IUPAC_equals(a1,a2)
    return IUPACnucLookup[a1,a2]
end

"""Needleman-Wunch, with IUPAC codes, and end gaps penalized slightly less."""
function IUPAC_nuc_edit_dist(s1::String, s2::String; edge_reduction = 0.9999)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 0.0


    s1numArr = toIUPACnum(s1)
    s2numArr = toIUPACnum(s2)
    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    arr = zeros(length(s1arr)+1, length(s2arr)+1)
    traceArr = zeros(Int, length(s1arr), length(s2arr))

    # this will need to be generalized when we want to allow overhang
    arr[:, 1] = edge_reduction*del_cost*(0:length(s1arr))
    arr[1,:] = edge_reduction*ins_cost*(0:length(s2arr))

    for i in 2:length(s1arr)+1
        for j in 2:length(s2arr)+1
            if s1arr[i-1] == s2arr[j-1]
                diag = arr[i-1, j-1] + match_cost
            elseif IUPAC_equals(s1numArr[i-1],s2numArr[j-1])
                diag = arr[i-1, j-1] + match_cost-0.00001
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = 1
            if i == length(s1arr)+1
                delMult = edge_reduction
            end
            insMult=1
            if j == length(s2arr)+1
                insMult = edge_reduction
            end

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            scores = [diag, del, ins]
            best = findmax(scores)[2]
            if best == 1
                arr[i, j] = diag
            elseif best ==2
                arr[i, j] = del
            elseif best ==3
                arr[i, j] = ins
            end
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]
    return alignedScore
end

"""Needleman-Wunch, with IUPAC codes, and end gaps penalized slightly less."""
function IUPAC_nuc_nw_align(s1::String, s2::String; edge_reduction = 0.9999)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 0.0


    s1numArr = toIUPACnum(s1)
    s2numArr = toIUPACnum(s2)
    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    arr = zeros(length(s1arr)+1, length(s2arr)+1)
    traceArr = zeros(Int, length(s1arr), length(s2arr))

    # this will need to be generalized when we want to allow overhang
    arr[:, 1] = edge_reduction*del_cost*(0:length(s1arr))
    arr[1,:] = edge_reduction*ins_cost*(0:length(s2arr))

    for i in 2:length(s1arr)+1
        for j in 2:length(s2arr)+1
            if s1arr[i-1] == s2arr[j-1]
                diag = arr[i-1, j-1] + match_cost
            elseif IUPAC_equals(s1numArr[i-1],s2numArr[j-1])
                diag = arr[i-1, j-1] + match_cost-0.00001
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = 1
            if i == length(s1arr)+1
                delMult = edge_reduction
            end
            insMult=1
            if j == length(s2arr)+1
                insMult = edge_reduction
            end

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            scores = [diag, del, ins]
            best = findmax(scores)[2]
            if best == 1
                arr[i, j] = diag
            elseif best ==2
                arr[i, j] = del
            elseif best ==3
                arr[i, j] = ins
            end
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]

    # return arr

    # the trace endpoing will need to be generalized if we want to
    # allow overhang.
    trI, trJ = length(s1arr), length(s2arr)

    # First compute the trace running backwards. Initialized to the
    # maximum size. With unit penalties, is it possible to tell how
    # big this should be in advance?
    backtrace = Array{Int}(undef, length(s1arr) + length(s2arr))
    # If you hit any boundary, you have to run all the way to the side.
    btInd = 1
    while (trI > 0) && (trJ > 0)
        backtrace[btInd] = traceArr[trI, trJ]
        if backtrace[btInd] == 1
            trI += -1
            trJ += -1
        elseif backtrace[btInd] == 2
            trJ += -1
        elseif backtrace[btInd] == 3
            trI += -1
        end
        btInd += 1
    end

    # If you hit the boundaries not at the top left corner.
    while trI > 0
        backtrace[btInd] = 3
        btInd += 1
        trI += -1
    end

    while trJ > 0
        backtrace[btInd] = 2
        btInd += 1
        trJ += -1
    end

    backtrace = backtrace[1:btInd-1]

    # This will need to be generalized to work on non-strings. Not
    # important.
    ali1arr = Array{Char}(undef, length(backtrace))
    ali2arr = Array{Char}(undef, length(backtrace))
    ind1 = 1
    ind2 = 1
    for i in 1:length(backtrace)
        tr = backtrace[length(backtrace)-(i-1)]
        if tr == 1
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = s2arr[ind2]
            ind1 += 1
            ind2 += 1
        elseif tr == 2
            ali2arr[i] = s2arr[ind2]
            ali1arr[i] = '-'
            ind2 += 1
        elseif tr == 3
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = '-'
            ind1 += 1
        end
    end
    return join(ali1arr), join(ali2arr), alignedScore
end


function banded_edit_dist(s1::String, s2::String; edge_reduction = 0.99, band_coeff = 1)
    # calculate band width
    avg_len = (length(s1) + length(s2)) / 2
    bandwidth = Int(ceil(band_coeff * sqrt(avg_len)))
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 1.0

    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal

    # dimension difference added to band length on one side
    # to give padding for large sequence length differences
    dim_diff = length(s1arr) - length(s2arr)
    maxlen = max(length(s1), length(s2))
    arr = zeros(maxlen + 1, 2*bandwidth + 1 + abs(dim_diff))
    traceArr = zeros(Int, size(arr))
    for i in 1:bandwidth+1
        add_to_band!(arr, edge_reduction*del_cost*(i-1), i, 1, bandwidth, dim_diff)
        add_to_band!(arr, edge_reduction*ins_cost*(i-1), 1, i, bandwidth, dim_diff)
    end

    for i in 2:length(s1arr)+1
        # keep j within 2*bandwidth + 1 band of i, but keep bandwidth sized buffer on edges
        lo = 0
        # to keep a big enough padding around the edges when there's a big length mismatch
        if i <= bandwidth + 1 + max(0, dim_diff)
            lo = 2
        else
            lo = i - bandwidth - max(0, dim_diff)
        end
        for j in lo:length(s2arr)+1
            # keep j within band
            if !(in_band(i, j, bandwidth, dim_diff))
                break
            end

            if in_band(i-1, j-1, bandwidth, dim_diff) && s1arr[i-1] == s2arr[j-1]
                diag = get_band_val(arr, i-1, j-1, bandwidth, dim_diff) + match_cost
            else
                diag = get_band_val(arr, i-1, j-1, bandwidth, dim_diff) + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = (i == length(s1arr)+1) ? edge_reduction : 1
            insMult = (j == length(s2arr)+1) ? edge_reduction : 1

            ins = in_band(i-1, j, bandwidth, dim_diff) ?
                    (get_band_val(arr, i-1, j, bandwidth, dim_diff)+(ins_cost*insMult)) : NaN
            del = in_band(i, j-1, bandwidth, dim_diff) ?
                    (get_band_val(arr, i, j-1, bandwidth, dim_diff)+(del_cost*delMult)) : NaN

            scores = [diag, del, ins]
            best = findmax(scores)[2]
            add_to_band!(arr, scores[best], i, j, bandwidth, dim_diff)
            add_to_band!(traceArr, best, i-1, j-1, bandwidth, dim_diff)
        end
    end
    alignedScore = get_band_val(arr, length(s1arr)+1, length(s2arr)+1, bandwidth, dim_diff)

    return alignedScore
end


#Indmax in 2D, returning the coordinates
function maxcoord(mat)
    return findmax(mat)[2]
end


#Dictionary-based matching
function ambig_positions(seq::String)
    return [1:length(seq);][[!(s in ['A','C','G','T']) for s in seq]]
end
function ambig_positions(seq::Array{Char})
    return ambig_positions(join(seq))
end
function IUPAC_char_equals(c1,c2)
    return IUPAC_equals(IUPACnumDict[c1],IUPACnumDict[c2])
end
function ambig_options(c::Char)
    return ['A','C','G','T'][[IUPAC_char_equals('A',c),IUPAC_char_equals('C',c),IUPAC_char_equals('G',c),IUPAC_char_equals('T',c)]]
end
function replace_positions(seq,pos,chars)
    chararr = collect(seq)
    for i in 1:length(pos)
        chararr[pos[i]] = chars[i]
    end
    return join(chararr)
end
function ambig_expand(seq)
    inds = ambig_positions(seq)
    options = [ambig_options(seq[i]) for i in inds]
    collec = [[]]
    for to_add in options
        collec = vcat([[vcat(el,[addit]) for addit in to_add] for el in collec]...)
    end
    return [replace_positions(seq,inds,coll) for coll in collec]
end
function all_mutants_of_all_ambigs(seq)
    ambig_seqs = ambig_expand(seq)
    return vcat([all_mutants(s) for s in ambig_seqs]...)
end
function mutate_arr_at_pos!(arr,pos,elem)
    arr[pos] = elem
    return arr
end
function all_mutants(seq)
    arr_seq = collect(seq)
    delsA = join.([vcat(deleteat!(copy(arr_seq),i),['A']) for i in 1:length(arr_seq)])
    delsC = join.([vcat(deleteat!(copy(arr_seq),i),['C']) for i in 1:length(arr_seq)])
    delsG = join.([vcat(deleteat!(copy(arr_seq),i),['G']) for i in 1:length(arr_seq)])
    delsT = join.([vcat(deleteat!(copy(arr_seq),i),['T']) for i in 1:length(arr_seq)])
    insertsA = join.([insert!(copy(arr_seq),i,'A')[1:end-1] for i in 1:length(arr_seq)])
    insertsC = join.([insert!(copy(arr_seq),i,'C')[1:end-1] for i in 1:length(arr_seq)])
    insertsG = join.([insert!(copy(arr_seq),i,'G')[1:end-1] for i in 1:length(arr_seq)])
    insertsT = join.([insert!(copy(arr_seq),i,'T')[1:end-1] for i in 1:length(arr_seq)])
    mutatesA = join.([mutate_arr_at_pos!(copy(arr_seq),i,'A')[1:end] for i in 1:length(arr_seq)])
    mutatesC = join.([mutate_arr_at_pos!(copy(arr_seq),i,'C')[1:end] for i in 1:length(arr_seq)])
    mutatesG = join.([mutate_arr_at_pos!(copy(arr_seq),i,'G')[1:end] for i in 1:length(arr_seq)])
    mutatesT = join.([mutate_arr_at_pos!(copy(arr_seq),i,'T')[1:end] for i in 1:length(arr_seq)])
    return vcat([seq],delsA,delsC,delsG,delsT,insertsA,insertsC,insertsG,insertsT,mutatesA,mutatesC,mutatesG,mutatesT)
end

function fast_primer_match(seqs,primers; tol_one_error = true)
    #Only tolerates a single bp difference between primer and seq
    #Assumes all filtering primers are the same length.
    #Note this doesn't mean that the primers actually had to be the same length!
    l = length(primers[1])
    matches = zeros(Int,length(seqs))
    noisy_primer_map = Dict{String,Int}()
    for (i,pr) in enumerate(primers)
        if tol_one_error
            targets = all_mutants_of_all_ambigs(pr)
        else
            targets = ambig_expand(pr)
        end
        for npr in targets
            noisy_primer_map[npr] = i
        end
    end
    for (i,s) in enumerate(seqs)
        if matches[i] == 0
            matches[i] = get(noisy_primer_map,s[1:l],0)
        end
        if matches[i] == 0
            matches[i] = -get(noisy_primer_map,reverse_complement(s[end-l+1:end]),0)
        end
    end
    return matches #Return sign is the rev_comp direction. Magnitude is the primer index.
end

function fast_primer_pair_match(seqs,fwd_primers,rev_primers; tol_one_error = true)
    fwd_matches = fast_primer_match(seqs,fwd_primers,tol_one_error=tol_one_error)
    rev_matches = fast_primer_match(seqs,rev_primers,tol_one_error=tol_one_error)
    problem_count = sum(fwd_matches.*rev_matches .> 0)
    if problem_count>0
        @warn "Inconsistencies: $(problem_count)"
    end
    keepers = fwd_matches.*rev_matches .< 0
    rev_comp_bool = fwd_matches .< 0
    return keepers,abs.(fwd_matches),abs.(rev_matches),rev_comp_bool
end

function demux_dict(seqs,fwd_primers,rev_primers; verbose = false, phreds = nothing, tol_one_error = true)
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
            end
        end
        return seq_dict
    end
end

function primer_trim(seq,primer; buffer = 3)
    a1,a2,score = IUPAC_nuc_nw_align(primer,seq[1:length(primer)+buffer], edge_reduction = 0.5)
    gapbool = reverse(collect(a1)) .!= '-'
    return seq[length(primer)-findfirst(gapbool)+buffer+2:end]
end
function primer_trim(seq,phreds,primer; buffer = 3)
    a1,a2,score = IUPAC_nuc_nw_align(primer,seq[1:length(primer)+buffer], edge_reduction = 0.5)
    gapbool = reverse(collect(a1)) .!= '-'
    return seq[length(primer)-findfirst(gapbool)+buffer+2:end],phreds[length(primer)-findfirst(gapbool)+buffer+2:end]
end

function primer_trim_reverse(seq,primer; buffer = 3)
    seq_tail = reverse_complement(seq[end-(length(primer)+buffer)+1:end])
    a1,a2,score = IUPAC_nuc_nw_align(primer,seq_tail, edge_reduction = 0.5)
    gapbool = reverse(collect(a1)) .!= '-'
    return seq[1: end-(length(primer)-findfirst(gapbool)+buffer+1)]
end
function primer_trim_reverse(seq,phreds,primer; buffer = 3)
    seq_tail = reverse_complement(seq[end-(length(primer)+buffer)+1:end])
    a1,a2,score = IUPAC_nuc_nw_align(primer,seq_tail, edge_reduction = 0.5)
    gapbool = reverse(collect(a1)) .!= '-'
    return seq[1: end-(length(primer)-findfirst(gapbool)+buffer+1)], phreds[1: end-(length(primer)-findfirst(gapbool)+buffer+1)]
end

function double_primer_trim(seq,fwd_primer,rev_primer; buffer = 3)
    t1 = primer_trim(seq,fwd_primer, buffer = buffer)
    return primer_trim_reverse(t1,rev_primer, buffer = buffer)
end
function double_primer_trim(seq,phreds,fwd_primer,rev_primer; buffer = 3)
    t1,p1 = primer_trim(seq,phreds,fwd_primer, buffer = buffer)
    return primer_trim_reverse(t1,p1,rev_primer, buffer = buffer)
end

function primer_peek(seqs::Vector{String}; l = 20, N = 30, keep = 20)
    #Assumptions: No ambigs, and primer sequence goes "barcode+primer"
    starts = [s[1:l] for s in seqs]
    ends = [s[end-l+1:end] for s in seqs]
    pot_primers = countmap(vcat(starts,reverse_complement.(ends)))
    topN = reverse(sort([(pot_primers[k],k) for k in keys(pot_primers)]))[1:N];
    for i in enumerate(topN)
        println(i)
    end
    return reverse_complement.(sort(reverse_complement.([i[2] for i in topN[1:keep]])))
end


#------Chunked demux funcions--------

"""
function chunked_fastq_filter_demux(chunk, chunk_size, seqs, phreds, names;
    demux_dir = "demux", fwd_primers = String[], rev_primers = String[],
    verbose = false, tol_one_error = true, primer_lookup = Dict(), error_rate = 0.01,
    min_length = 30, max_length = 1000000, label_prefix = "seq", error_out = true)

Intended to be passed to `chunked_fastq_apply()` for quality/length filtering and demultiplexing of
FASTQ files in one pass.
"""
function chunked_fastq_filter_demux(chunk, chunk_size, seqs, phreds, names;
    demux_dir = "demux", rejects_dir = "rejects", fwd_primers = String[], rev_primers = String[],
    verbose = false, tol_one_error = true, primer_lookup = Dict(), error_rate = 0.01,
    min_length = 30, max_length = 1000000, label_prefix = "seq", error_out = true)
    
    total_reads, quality_reads, demuxed_reads = 0, 0, 0
    total_reads = length(seqs)
    
    if verbose
        println("processing chunk $(chunk) of size $(length(seqs))")
    end

    #filter...
    lengths = length.(seqs)
    mean_errors = [mean(phred_to_p.(phred)) for phred in phreds]
    inds = [1:length(seqs);][(lengths .< max_length) .& (lengths .> min_length) .& (mean_errors .< error_rate)]
    rejects = setdiff(1:length(seqs),inds)
                                                                                            
    # write failed sequences
    outpath = rejects_dir*"/QUALITY_REJECTS.fastq"
    records = FASTQ.Record.(
                    names[rejects],
                    seqs[rejects],
                    phreds[rejects]
                    )
    writer = FASTQ.Writer(open(outpath, "a"))
    for record in records
         write(writer, record)
    end
    close(writer)
                                                                                            
    if error_out == true
        filtered_names = ["$label_prefix$(((chunk-1)*chunk_size+i))|ee=$(mean_errors[i])" for i in inds]
    else
        filtered_names = ["$label_prefix$(((chunk-1)*chunk_size+i))" for i in 1:length(inds)]
    end

    filtered_seqs, filtered_phreds = seqs[inds], phreds[inds]
    quality_reads = length(filtered_seqs)

    #demux...
    demux_dic = demux_dict(filtered_seqs, fwd_primers, rev_primers;
        phreds = filtered_phreds, tol_one_error = tol_one_error
    )
    #write...
    for (i, fwd_primer) in enumerate(fwd_primers)
        for (j, rev_primer) in enumerate(rev_primers)
            if haskey(demux_dic,(i,j))
                if haskey(primer_lookup,(fwd_primer,rev_primer))
                    ID = primer_lookup[fwd_primer,rev_primer]
                    outpath = demux_dir*"/$(ID)_$(fwd_primer)_$(rev_primer).fastq"
                else
                    ID = "UNKNOWN"
                    outpath = rejects_dir*"/$(ID)_$(fwd_primer)_$(rev_primer).fastq"
                end
                result = demux_dic[i,j]
                records = FASTQ.Record.(
                    filtered_names[getindex.(result, 3)],
                    getindex.(result, 1),
                    getindex.(result, 2)
                    )
                writer = FASTQ.Writer(open(outpath, "a"))
                for record in records
                    write(writer, record)
                    if ID != "UNKNOWN"
                        demuxed_reads += 1
                    end
                end
                close(writer)
            end
        end
    end
    return [total_reads, quality_reads, demuxed_reads]
end
