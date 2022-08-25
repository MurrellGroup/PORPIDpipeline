using NextGenSeqUtils, PORPID, StatsBase, HypothesisTests, DataFrames,
    BioSequences, IterTools, CSV, FASTX

"""
Given results of the LDA, filter likely offspring UMIs and apply qc
"""
function filterCCSFamilies(most_likely_real_for_each_obs, path, index_to_tag, tag_counts, template_name, template; fs_thresh=5, lda_thresh=0.995)

    # first check if there are any REJECTS and if there are annotate them with
    # the tag BPB-rejects and copy them to reject file one level up
    run(`touch $(path)REJECTS.fastq`)
    seqs, phreds, names = read_fastq(path*"REJECTS.fastq")
    write_fastq(path[1:end-1]*"_rejects.fastq", seqs, phreds,
                     names=names .* " BPB-rejects", append=true)
    run(`rm $(path)REJECTS.fastq`)
    
    likely_real = []
    UMI = []
    bin_sizes = []
    tags = []
    probs = []

    for (observed_index_local, tuple_local) in enumerate(most_likely_real_for_each_obs)
        tagSeq = index_to_tag[observed_index_local]
        prob = tuple_local[2]

        push!(UMI, tagSeq)
        push!(bin_sizes, tag_counts[tagSeq])
        push!(probs, log(1-prob))

        #Testing for whether there is a drop in QC scores in the barcode region
        #This filters out sequencing errors, heteroduplexes
        read_seqs, read_phreds, read_names = read_fastq(path*tagSeq*".fastq");
        averages = mean([Array{Float64}(read_phreds[i][1:50]) for i in 1:length(read_phreds)])
        umi_ix = findfirst(r"n+", template)
        p_val = 1
        if mean(averages[umi_ix])<mean(averages[umi_ix[end]+5:umi_ix[end]+25])
            p_val=pvalue(EqualVarianceTTest(averages[umi_ix],averages[umi_ix[end]+5:umi_ix[end]+25]))
        end

        if p_val < 1/(5*tag_counts[tagSeq]^2) && minimum(averages[umi_ix]) < mean(averages[umi_ix[end]+5:umi_ix[end]+25])/2
            push!(tags, "heteroduplex")
            write_fastq(path[1:end-1]*"_rejects.fastq", read_seqs, read_phreds,
                      names=read_names .* " heterodupex", append=true)
            continue
        end

        if prob < lda_thresh # used to be .995
            push!(tags, "LDA-rejects")
            write_fastq(path[1:end-1]*"_rejects.fastq", read_seqs, read_phreds,
                    names=read_names .* " LDA-rejects", append=true)
            continue
        end

        #Filter out if not length 8
        if length(tagSeq) != 8
            push!(tags, "UMI_len != 8")
            write_fastq(path[1:end-1]*"_rejects.fastq", read_seqs, read_phreds,
                    names=read_names .* "UMI_len != 8", append=true)
            continue
        end

        # fs_thresh =  # used to be 5
        if tag_counts[tagSeq] < fs_thresh # set fs_thresh here
            push!(tags, "fs<$(fs_thresh)")
            write_fastq(path[1:end-1]*"_rejects.fastq", read_seqs, read_phreds,
                    names=read_names .* "fs<$(fs_thresh)", append=true)
            continue
        end

        push!(tags, "likely_real")
        push!(likely_real, (prob, tagSeq))
    end

    #sort!(likely_real, by=x->tag_counts[x[2]])
        directory = "$(path[1:end-1])_keeping/"
    mkdir(directory)

    for tag_tuple in likely_real
        run(`cp $(path)$(tag_tuple[2]).fastq $(directory)$(tag_tuple[2])_$(length(tag_tuple[2]))_$(round(tag_tuple[1],digits = 5))_$(tag_counts[tag_tuple[2]]).fastq`)
    end
    tag_df = DataFrame(Sample=template_name, UMI = UMI, fs = bin_sizes, tags = tags, probs = probs)
    return tag_df
end

#writing PORPID output
function porpid_write_to_file(source_file_name, template, tag, output_sequence,
  score, outdir)
    source_file_name = basename(source_file_name)
    output_file_name = "$(outdir)/$(source_file_name)/$(template.name)/$(tag).fastq"
    mkpath("$(outdir)/$(source_file_name)/$(template.name)")
    fo = open(output_file_name, "a")
    writer = FASTQ.Writer(fo)
    write(writer, output_sequence)
    close(writer)
    close(fo)
end

function porpid_write_to_dictionary(dictionary, source_file_name, template, tag,
  output_sequence, score)
    directory = "$(source_file_name)/$(template.name)"
    if !haskey(dictionary, directory)
        dictionary[directory] = Dict()
    end
    directory_dict = dictionary[directory]
    if !haskey(directory_dict, tag)
        directory_dict[tag] = []
    end
    push!(directory_dict[tag], (score, output_sequence))
end

function porpid_write_to_file_count_to_dict(dictionary, source_file_name, template,
  tag, output_sequence, score, outdir)
    porpid_write_to_file(source_file_name, template, tag, output_sequence, score, outdir)
    directory = "$(source_file_name)/$(template.name)"
    if !haskey(dictionary, directory)
        dictionary[directory] = Dict()
    end
    directory_dict = dictionary[directory]
    directory_dict[tag] = get(directory_dict, tag, 0) + 1
end

