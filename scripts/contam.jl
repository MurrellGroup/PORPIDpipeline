ENV["MPLBACKEND"] = "Agg"
using PorpidPostproc, NextGenSeqUtils, BioSequences, DataFrames,CSV, DataFramesMeta
# include("../../src/functions.jl")
# include("../../src/contam-filter_functions.jl")

proportion_thresh = snakemake.params["proportion_thresh"]
cluster_thresh = snakemake.params["cluster_thresh"]
dist_thresh = snakemake.params["dist_thresh"]
contam_toggle = snakemake.params["contam_toggle"]
files = snakemake.input["files"]


run_ID = snakemake.wildcards["dataset"]
contam_passed_dir = snakemake.output[1]
contam_failed_dir = snakemake.output[2]
mkpath(contam_passed_dir)
mkpath(contam_failed_dir)

contam_df = DataFrame(  sample = String[],
                        sequence_name = String[],
                        nearest_nonself_variant = String[],
                        nearest_nonself_distance = Float64[],
                        flagged = Bool[],
                        discarded = Bool[]);
                        
suspect_contam_df = DataFrame(  sample = String[], sequence_name = String[],
            nearest_nonself_variant = String[], nearest_nonself_distance = Float64[],
            flagged = Bool[]);
                        
if contam_toggle == "on"
println("performing contam filter...")

    one_run_db = vcat([db_cluster_seqs(f,
                            proportion_thresh = proportion_thresh,
                            cluster_thresh = cluster_thresh) for f in files]...);
    contam_db = db_seqs(snakemake.input["panel"])
    merged_db = vcat(one_run_db,contam_db);

    for f in files
        sample = split(basename(f),".fast")[1] #modify this if you want to get the dataset without the file extension.
        flagged_names,flagged_neighbours,flagged_dists,reported,discarded,discarded_bool_inds = contam_check(f,merged_db, thresh=dist_thresh)
        seqnames,seqs = read_fasta_with_descriptors_in_names(f)
        if length(flagged_names) > 0
            for r in 1:length(flagged_names)
                #println((flagged_names,flagged_neighbours,flagged_dists,reported,discarded))
            push!(contam_df,[sample,flagged_names[r],basename(flagged_neighbours[r]),flagged_dists[r],reported[r],discarded[r]])
            end
            write_fasta(contam_passed_dir*"/"*sample*".fasta",seqs[.!discarded_bool_inds],names = seqnames[.!discarded_bool_inds])
            write_fasta(contam_failed_dir*"/"*sample*".fasta",seqs[discarded_bool_inds],names = seqnames[discarded_bool_inds])
        else
            write_fasta(contam_passed_dir*"/"*sample*".fasta",seqs,names = seqnames)
        end
    end
    CSV.write(snakemake.output[3],contam_df)

    # do this again with proportion_thresh = 0.0 to find further contam suspects
    suspect_run_db = vcat([db_cluster_seqs(f,
                        proportion_thresh = 0.0,
                        cluster_thresh = cluster_thresh) for f in files]...);
    suspect_merged_db = vcat(suspect_run_db,contam_db);
    
    for f in files
        sample = split(basename(f),".fast")[1]
        flagged_names, flagged_neighbours, flagged_dists, reported, suspects, suspects_bool_inds = contam_check(f,suspect_merged_db, thresh=dist_thresh)
        seqnames,seqs = read_fasta_with_descriptors_in_names(f)
        if length(flagged_names) > 0
            for r in 1:length(flagged_names)
                if ( reported[r] )
                    push!(suspect_contam_df, [sample, flagged_names[r], basename(flagged_neighbours[r]), flagged_dists[r], reported[r]])
                end
            end
        end
    end
    sort!(suspect_contam_df,[:nearest_nonself_distance])
    CSV.write(snakemake.output[4],suspect_contam_df)
else
    println("bypassing contam filter...")
    for f in files
        sample = split(basename(f),".fast")[1]
        seqnames,seqs = read_fasta_with_descriptors_in_names(f)
        write_fasta(contam_passed_dir*"/"*sample*".fasta",seqs,names = seqnames)
    end
    CSV.write(snakemake.output[3],contam_df)
    CSV.write(snakemake.output[4],suspect_contam_df)
end
