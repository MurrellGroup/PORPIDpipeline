using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, RobustAmpliconDenoising, NextGenSeqUtils, CSV, DataFrames

function generateConsensusFromDir(dir, template_name)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        exit()
    end
    cons_collection = map(ConsensusFromFastq, files)
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    return seq_collection, seqname_collection
end

function ConsensusFromFastq(file)
    seqs,phreds,seq_names = read_fastq(file)
    draft = consensus_seq(seqs)
    draft2 = refine_ref(draft, seqs)
    final_cons = refine_ref(draft2,seqs)
    alignments, maps, matches, matchContent = getReadMatches(final_cons, seqs, 0)
    cons_name = split(basename(file),"_")[1]*" num_CCS=$(length(seqs)) min_agreement=$(round(minimum(matches); digits = 2))"
    return final_cons, cons_name
end

"""
Returns an array of degapped coordinates, such that
coords(ref, read)[i] gives you the position the aligned read/ref
that matches the i'th ungapped position in ref.
"""
function coords(ref, read)
    if length(ref) != length(read)
        error("Aligned strings are meant to be the same length.")
    end
    degappedRef = degap(ref)
    coordMap = zeros(Int64, length(degappedRef))
    count = 1
    for i in 1:length(degappedRef)
        while ref[count] == '-'
            count += 1
        end
        coordMap[i] = count
        count += 1
    end
    return coordMap
end

"""
Return matches to a candidate reference from a set of reads.
"""
function getReadMatches(candidate_ref, reads, shift; degap_param = true, kmer_align = true)
    alignments = []
    if kmer_align
        alignments = map(i -> kmer_seeded_align(candidate_ref, i), reads)
    else
        alignments = map(i -> nw_align(candidate_ref, i), reads)
    end

    maps = [coords(i...) for i in alignments]

    if (degap_param)
        matchContent = [[degap(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
        matches = [freq(matchContent[k], degap(candidate_ref[k:k+shift])) for k in 1:length(matchContent)]
    else
       matchContent = [[(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
       matches = [freq(matchContent[k], candidate_ref[k:k+shift]) for k in 1:length(matchContent)]
    end
    return alignments, maps, matches, matchContent
end

config = snakemake.params["config"]
#Calculate consensus sequences for each family.
t1 = time()
template_name = snakemake.wildcards["sample"]
cDNA_primer = config["cDNA_primer"]
SID_ix = findfirst(r"[a-z]+", cDNA_primer)
to_trim = uppercase(cDNA_primer[SID_ix[1]:end])
println("Processing $(template_name)")
base_dir = snakemake.input[1]*"/"*template_name*"_keeping"
seq_collection, seqname_collection = generateConsensusFromDir(base_dir, template_name)
trimmed_collection = [primer_trim(s,to_trim) for s in seq_collection];
write_fasta(snakemake.output[1],reverse_complement.(trimmed_collection),names = seqname_collection)

# Update tag data with minimum_agreement and artefact rejects
tag_df = CSV.read(snakemake.input[2], DataFrame)

# first do min agrement filter
agreement_thresh = snakemake.params["agreement_thresh"]
minagrs=(x->parse(Float64,split(split(x," ")[3],"=")[2])).(seqname_collection)
minagr_rejects = (x->split(x," ")[1][end-7:end]).(seqname_collection[minagrs .< agreement_thresh])
minag_count=0
for row in eachrow(tag_df)
    if row[:tags] == "likely_real" && row[:UMI] in minagr_rejects
        row[:tags]="minag-reject"
        global minag_count+=1
    end
end
println("$(template_name): labelling $(minag_count) reads as minag-reject")


# now rename possible artefacts
af_thresh = snakemake.params["af_thresh"]
ccs=tag_df[tag_df[!,:tags].=="likely_real",:fs]
af_cutoff=artefact_cutoff(ccs,af_thresh)
art_count=0
for row in eachrow(tag_df)
   if row[:tags] == "likely_real" && row[:fs]<af_cutoff
        row[:tags]="maybe-artefact"
        global art_count+=1
   end
end
println("$(template_name): labeling $(art_count) reads with fs under $(af_cutoff) as maybe-artefact")

CSV.write(snakemake.output[2], sort!(tag_df, [:Sample, :tags, :fs], rev = [false, false, true]));

t2 = time()
println("Consensus generation took $(t2-t1) seconds.")
