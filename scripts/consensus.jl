using Distributed
@everywhere ENV["MPLBACKEND"] = "Agg"
@everywhere using RobustAmpliconDenoising, NextGenSeqUtils

function generateConsensusFromDir(dir, template_name)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq"]
    if length(files) > 0
        println("Generating consensus for $(length(files)) templates")
    else
        println("WARNING: no template families for $(template_name)")
        exit()
    end
    cons_collection = pmap(ConsensusFromFastq, files)
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    return seq_collection, seqname_collection
end

@everywhere function ConsensusFromFastq(file)
    seqs,phreds,seq_names = read_fastq(file)
    draft = consensus_seq(seqs)
    draft2 = refine_ref(draft, seqs)
    final_cons = refine_ref(draft2,seqs)
    alignments, maps, matches, matchContent = getReadMatches(final_cons, seqs, 0)
    cons_name = split(basename(file),"_")[1]*" num_CCS=$(length(seqs)) min_agreement=$(round(minimum(matches); digits = 2))"
    return final_cons, cons_name
end

@everywhere begin
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
end

@everywhere begin
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
t2 = time()
println("Consensus generation took $(t2-t1) seconds.")
