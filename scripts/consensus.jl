using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, CSV, DataFrames, StatsBase, BioSequences

import RobustAmpliconDenoising

ag_directory = snakemake.output[3]

function generateConsensusFromDir(dir, template_name)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating $(template_name) consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        exit()
    end
    cons_collection = map(ConsensusFromFastq, files)
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    return seq_collection, seqname_collection
end

function runLengths(vec)
    (vals,lens)=rle(vec)
    ret=vcat([ones(Int,lens[i])*lens[i] for i in 1:length(lens)]...)
    return ret
end

function complement_char(s::String)
    c='-'
    length(s)>0 ? c=s[1] : nothing
    return string( complement(DNA(c)) )
end

function ConsensusFromFastq(file)
    ma_margin=1 # to ignore first and last nt in min_ag calculation.
    seqs,phreds,seq_names = PORPIDpipeline.read_fastq(file)
    draft = RobustAmpliconDenoising.consensus_seq(seqs)
    draft2 = RobustAmpliconDenoising.refine_ref(draft, seqs)
    final_cons = RobustAmpliconDenoising.refine_ref(draft2,seqs)
    alignments, maps, matches, matchContent = getReadMatches(final_cons, seqs, 0)
    min_ag=round(minimum(matches[1+ma_margin:end-ma_margin]); digits = 2)
    cons_name = split(basename(file),"_")[1]*" fs=$(length(seqs)) minag=$(min_ag)"
    ag_list = [round(m; digits = 2) for m in matches] #make list of agreement scores
    # min_ag = minimum(ag_list)
    ag_nts = (complement_char).( (mode).(matchContent) )
    ag_rls = runLengths(ag_nts)
    ag_inds = (x->x<=min_ag).(ag_list)
    ag_inds = ag_inds .& (x->x!="-").(ag_nts) # drop gap characters from min ag dataframe
    ag_nt = ag_nts[ag_inds]
    ag_rl = ag_rls[ag_inds]
    ag_pos = (x->length(ag_list)-x+1).(collect(1:length(ag_list))[ag_inds])
    ag_df = DataFrame(ag = ag_list[ag_inds], agp=ag_pos, agnt=ag_nt, agrl=ag_rl) #create dataframe from this data
    mkpath(ag_directory) #couldn't figure out why snakemake wouldn't make the path, had to do it here
    CSV.write(ag_directory*"/"*template_name*split(basename(file),"_")[1]*"_agreement.csv", ag_df);
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
        alignments = map(i -> RobustAmpliconDenoising.kmer_seeded_align(candidate_ref, i), reads)
    else
        alignments = map(i -> RobustAmpliconDenoising.nw_align(candidate_ref, i), reads)
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
trimmed_collection = [RobustAmpliconDenoising.primer_trim(s,to_trim) for s in seq_collection];
PORPIDpipeline.write_fasta(snakemake.output[1],PORPIDpipeline.reverse_complement.(trimmed_collection),
        names = seqname_collection)

# First update tag data with minimum_agreement scores
tag_df = CSV.read(snakemake.input[2], DataFrame)
tag_df.minag .= 0.0

minagrs=(x->( split(x," ")[1][end-7:end], parse(Float64,split(split(x," ")[3],"=")[2]) )).(seqname_collection)

# minagr_records = (x->split(x," ")[1][end-7:end]).(seqname_collection)
for row in eachrow(tag_df)
    if row[:tags] == "likely_real"
        for (umi,minag) in minagrs
            if row[:UMI] == umi
                row[:minag] = minag
            end
        end
    end
end
println("$(template_name): tag file updated with minag scores")


CSV.write(snakemake.output[2], sort!(tag_df, [:Sample, :tags, :fs], rev = [false, false, true]));

t2 = time()
println("Consensus generation for $(template_name) took $(t2-t1) seconds.")
