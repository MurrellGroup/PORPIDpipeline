using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, PORPID, StatsBase
using HypothesisTests, DataFrames, BioSequences, IterTools, CSV, FASTX


#iterate through samples, run PORPID, and filter families.
t1 = time()
config = snakemake.params["config"]
fs_thresh = snakemake.params["fs_thresh"]

# allow for an fs override
if "fs_override" in keys(config)
    fs_thresh = config["fs_override"]
end

lda_thresh = snakemake.params["lda_thresh"]
data_dir = snakemake.input[1]
ID = uppercase(match(r"[a-z]+", config["cDNA_primer"]).match)
umi_ix = findfirst(r"N+", config["cDNA_primer"])
#trying template with 6bp sample IDs and last 8bp primer
template_suffix = replace(config["cDNA_primer"][umi_ix[1]:end],
    "N" => "n")*"*"

filtered_data_file = snakemake.wildcards["sample"]*".fastq.gz"
templates = Dict()
templates[snakemake.wildcards["sample"]] = ID*template_suffix #fix

cfg = Configuration()
cfg.files = ["$(data_dir)/$(filtered_data_file)"]
cfg.filetype = fastq
cfg.start_inclusive = 0
cfg.end_inclusive = 6 + length(template_suffix) + 2 #setting this explicitly
cfg.try_reverse_complement = false #The sequences are already oriented
for (name, template) in templates
    push!(cfg.templates, Template(name, template))
end

println("$(snakemake.wildcards["sample"]): processing output/$(filtered_data_file)/...")
println("$(snakemake.wildcards["sample"]): using template $(ID*template_suffix)")
if isdir(snakemake.output[1])
    @error "The directory already exists! Please delete"
    #rm(snakemake.output[1], recursive=true)
end
dir_dict = Dict()
my_output_func(source_file_name,
    template,
    tag,
    output_sequence,
    score
    ) = porpid_write_to_file_count_to_dict(dir_dict,
                                           source_file_name,
                                           template,
                                           tag,
                                           output_sequence,
                                           score,
                                           dirname(snakemake.output[1])) #fix
say_print_func = function(count)
    println("$(snakemake.wildcards["sample"]): processed $(count) sequences")
end
# This is the slow bit
extract_tags_from_file(cfg.files[1],
    cfg,
    my_output_func,
    print_every=10000,
    print_callback=say_print_func)
directories = collect(keys(dir_dict))

#There is now only one template per demuxed dataset
analysing_template = 1;
template_name = basename(directories[analysing_template])
tag_dict = dir_dict[directories[analysing_template]]
BPB_rejects = 0
if "REJECTS" in keys(tag_dict)
    BPB_rejects = tag_dict["REJECTS"]
    # println(template_name," primer block REJECTS => ",tag_dict["REJECTS"])
    # parent_path = snakemake.output[1]*"/"
    # println(parent_path)
    # run(`mv $(parent_path)$(template_name)/REJECTS.fastq $(parent_path)REJECTS.fastq`)
end
println(template_name," BPB-rejects => ",BPB_rejects)
delete!(tag_dict, "REJECTS")
tag_counts = tag_dict
tags = collect(keys(tag_dict))

#Builts matrix of conditional probabilities.
tag_to_index, index_to_tag = tag_index_mapping(tags)
pacbio_error_rate = 0.005
recurse = 1
probabilities_array = prob_observed_tags_given_reals(tag_to_index,
                                                     PORPID.PacBioErrorModel(pacbio_error_rate),
                                                     recurse)
indexed_counts = index_counts(tag_counts, tag_to_index);
#Runs the LDA inference
most_likely_real_for_each_obs = LDA(probabilities_array, indexed_counts)
path = snakemake.output[1]*"/"*template_name*"/" #fix

#Filter and copy to "_keeping"
tag_df = filterCCSFamilies(most_likely_real_for_each_obs, path,
    index_to_tag,
    tag_counts,
    template_name,
    templates[template_name],
    fs_thresh=fs_thresh,
    lda_thresh=lda_thresh
)

if BPB_rejects > 0
    push!(tag_df, [template_name, "REJECTED", BPB_rejects, "BPB-rejects", -5.123456789])
end

CSV.write(snakemake.output[2], sort!(tag_df, [:Sample, :tags, :fs], rev = [false, false, true])); #fix
t2 = time()
println("UMI identification took $(t2-t1) seconds.")
