using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, CSV, BioSequences
using DataFrames, DataFramesMeta, StatsBase
import MolecularEvolution

fasta_collection = snakemake.input[1]*"/"*snakemake.wildcards["sample"]*".fasta"

sample = snakemake.wildcards["sample"]
dataset = snakemake.wildcards["dataset"]
fs_thresh = snakemake.params["fs_thresh"]
af_thresh = snakemake.params["af_thresh"]
q_thresh = snakemake.params["q_thresh"]
agreement_thresh = snakemake.params["agreement_thresh"]
panel_thresh = snakemake.params["panel_thresh"]

config = snakemake.params["config"]
if "fs_override" in keys(config)
    fs_thresh = config["fs_override"]
end

if "af_override" in keys(config)
    af_thresh = config["af_override"]
end

if "q_override" in keys(config)
    q_thresh = config["q_override"]
end

if "ma_override" in keys(config)
    agreement_thresh = config["ma_override"]
end

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


# update tag data with minimum_agreement and artefact rejects
tag_df = CSV.read(snakemake.input[2], DataFrame)

# first do min agrement filter
minag_count=0
for row in eachrow(tag_df)
    if row[:tags] == "likely_real" && row[:minag] < agreement_thresh
        row[:tags]="minag-reject"
        global minag_count+=1
    end
end
println("$(sample): labelling $(minag_count) families as minag-reject")


# now rename possible artefacts
ccs=tag_df[tag_df[!,:tags].=="likely_real",:fs]
af_cutoff=artefact_cutoff(ccs,af_thresh,q_thresh)
q_cutoff=Int(ceil(quantile(ccs,q_thresh))) # maximum_of_non_outliers(ccs,q_thresh)
art_count=0
for row in eachrow(tag_df)
   if row[:tags] == "likely_real" && row[:fs]<af_cutoff
        row[:tags]="maybe-artefact"
        global art_count+=1
   end
end
println("$(sample): labelling $(art_count) families with fs under $(af_cutoff) as maybe-artefact")

CSV.write(snakemake.output[11], sort!(tag_df, [:Sample, :tags, :fs], rev = [false, false, true]));

ali_seqs,seqnames,af_cutoff = H704_init_template_proc(fasta_collection, panel_file, snakemake.output[1], snakemake.output[2],  snakemake.output[3], snakemake.output[4],  agreement_thresh=agreement_thresh, panel_thresh=panel_thresh, af_thresh=af_thresh,q_thresh=q_thresh)

sp_selected = @linq tag_df |> where(:Sample .== sample)
sp_selected = @linq sp_selected |> where(:tags .!= "BPB-rejects")
fig = family_size_umi_len_stripplot(sp_selected,fs_thresh=fs_thresh,
        af_thresh=af_thresh, q_thresh=q_thresh,
        af_cutoff=af_cutoff, q_cutoff=q_cutoff)
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
        af_thresh=af_thresh, q_thresh=q_thresh,
        af_cutoff=af_cutoff, q_cutoff=q_cutoff)
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
  
sample_dir = snakemake.input[4]
sp_selected = @linq tag_df |> where(:Sample .== sample)
sp_minag_rejects = @linq sp_selected |> where(:tags .== "minag-reject")
sp_reals = @linq sp_selected |> where(:tags .== "likely_real")
sp_selected = vcat(sp_reals, sp_minag_rejects)
fig = minag_position_plot(sample_dir,sp_selected,ma_thresh=agreement_thresh)
fig.savefig(snakemake.output[12]; transparent = true, dpi = 200, bbox_inches = "tight")
