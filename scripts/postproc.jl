using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, CSV, NextGenSeqUtils, BioSequences, DataFrames, DataFramesMeta
import MolecularEvolution

# include("../../src/functions.jl")
# include("../../src/apobec_model.jl")
# include("../../src/molev_functions.jl")
# include("../../src/postproc_functions.jl")

fasta_collection = snakemake.input[1]*"/"*snakemake.wildcards["sample"]*".fasta"
tag_df = CSV.read(snakemake.input[2], DataFrame)
sample = snakemake.wildcards["sample"]
dataset = snakemake.wildcards["dataset"]
fs_thresh = snakemake.params["fs_thresh"]
af_thresh = snakemake.params["af_thresh"]
agreement_thresh = snakemake.params["agreement_thresh"]
panel_thresh = snakemake.params["panel_thresh"]

config = snakemake.params["config"]
if "fs_override" in keys(config)
    fs_thresh = config["fs_override"]
end

if "af_override" in keys(config)
    af_thresh = config["af_override"]
end

#env_seqs = read_fasta("panels/env_column_stripped_panel.fasta")
#env_profile = seqs2profile(uppercase.(env_seqs))

#re_seqs = read_fasta("panels/HIV1_COM_2017_5970-8795_DNA_stripped.fasta")
#re_profile = seqs2profile(uppercase.(re_seqs))

# check for fasta collection
if !isfile(fasta_collection)
    @error "No input FASTA file at $(fasta_collection)! Check if no output sequences for previous step."
    exit()
end

# check for panel file
if !isfile(snakemake.params["panel"])
    @error "Panel $(snakemake.params["panel"]) not found!"
    exit()
end
panel_file = snakemake.params["panel"]

# get af_cutoff from tags dataframe
# sp_selected = @linq tag_df |> where(:Sample .== sample)
# sp_artefacts = @linq sp_selected |> where(:tags .== "maybe-artefact")
# sp_reals = @linq sp_selected |> where(:tags .== "likely_real")
# sp_selected = vcat(sp_artefacts, sp_reals)
# fss = sp_selected[!,:fs]
# af_cutoff=artefact_cutoff(fss,af_thresh)

# fss = sp_selected[!,:fs]
# af_cutoff=1
# if length(fss)>0
#     af_cutoff=maximum(fss)+1
# end

ali_seqs,seqnames,af_cutoff = H704_init_template_proc(fasta_collection, panel_file, snakemake.output[1], snakemake.output[2],  snakemake.output[3], snakemake.output[4],  agreement_thresh=agreement_thresh, panel_thresh=panel_thresh, af_thresh=af_thresh)

sp_selected = @linq tag_df |> where(:Sample .== sample)
sp_selected = @linq sp_selected |> where(:tags .!= "BPB-rejects")
fig = family_size_umi_len_stripplot(sp_selected,fs_thresh=fs_thresh,
        af_thresh=af_thresh,af_cutoff=af_cutoff)
fig.savefig(snakemake.output[5];
    transparent = true,
    dpi = 200,
    bbox_inches = "tight")
    
sp_selected = @linq tag_df |> where(:Sample .== sample)
sp_artefacts = @linq sp_selected |> where(:tags .== "maybe-artefact")
sp_minag_rejects = @linq sp_selected |> where(:tags .== "minag-reject")
sp_fs_rejects = @linq sp_selected |> where(:tags .== "fs<$(fs_thresh)")
sp_reals = @linq sp_selected |> where(:tags .== "likely_real")
sp_selected = vcat(sp_artefacts, sp_reals, sp_minag_rejects, sp_fs_rejects)
fig = family_size_stripplot(sp_selected,fs_thresh=fs_thresh,
        af_thresh=af_thresh,af_cutoff=af_cutoff)
fig.savefig(snakemake.output[6];
    transparent = true,
    dpi = 200,
    bbox_inches = "tight")

selected = @linq tag_df |> where(:Sample .== sample)
gdf = DataFramesMeta.groupby(selected, :tags)
summary = @combine gdf cols(AsTable) = ( porpid_result=first(:tags), n_UMI_families=length(:fs), n_CCS=sum(:fs) )
# println( summary[!, [:porpid_result,:n_UMI_families,:n_CCS]] )
CSV.write(snakemake.output[7],summary[!, [:porpid_result,:n_UMI_families,:n_CCS]])

umi_dir=snakemake.input[3]*"/"*sample
umis = readdir(umi_dir)
umis = umis[(x -> findall(".",x)[1][1] > 0).(umis)]
weights = ones(length(umis))
for k in 1:length(umis)
  seqs, phreds, names = read_fastq(umi_dir*"/"*umis[k])
  weights[k] = length(seqs)
  j = findall(".",umis[k])[1][1]
  umis[k] = umis[k][1:j-1]
end
fig = di_nuc_freqs(umis, weights=weights )
fig.savefig(snakemake.output[10];
  transparent = true,
  dpi = 200,
  bbox_inches = "tight")
