using BioSequences, FASTX

const IUPACbool = Dict{Char,Array{Bool,1}}(Dict())
IUPACbool['A']=[true,false,false,false]
IUPACbool['C']=[false,true,false,false]
IUPACbool['G']=[false,false,true,false]
IUPACbool['T']=[false,false,false,true]
IUPACbool['U']=[false,false,false,true]
IUPACbool['R']=[true,false,true,false]
IUPACbool['Y']=[false,true,false,true]
IUPACbool['S']=[false,true,true,false]
IUPACbool['W']=[true,false,false,true]
IUPACbool['K']=[false,false,true,true]
IUPACbool['M']=[true,true,false,false]
IUPACbool['B']=[false,true,true,true]
IUPACbool['D']=[true,false,true,true]
IUPACbool['H']=[true,true,false,true]
IUPACbool['V']=[true,true,true,false]
IUPACbool['N']=[true,true,true,true];

function resolve_base(c)
    uc = uppercase(c)
    if uc in keys(IUPACbool)
        return ['A','C','G','T'][rand((1:4)[IUPACbool[uc]])]
    else
        return uc
    end
end

function resolve_seq(s)
    return join(resolve_base.(collect(s)))
end

function db_cluster_seqs(path; proportion_thresh = 0.2, cluster_thresh = 0.015)
    seqnames,seqs = read_fasta_with_names(path)
    seqs = degap.(seqs)
    resolved_seqs = resolve_seq.(seqs)
    kmer_vecs = kmer_count.(resolved_seqs,6)
    μs, sizes, cluster_indices, centroid_indices = dp_centers(kmer_vecs, cluster_thresh,
                                                              distfunc=corrected_kmer_dist(6),
                                                              center=mean)
    print(".")
    props = sizes ./ sum(sizes)
    keeps = props .> proportion_thresh
    for_return = collect(zip([path for i in keeps[keeps]].*"_".*string.(round.(100 .* props[keeps],sigdigits = 3)).*"%",μs[keeps],seqs[centroid_indices[keeps]]))
    push!(for_return,(path*"_All",mean(kmer_vecs),mode(seqs)))
    return for_return
end

function db_seqs(path)
    seqnames,seqs = read_fasta_with_names(path)
    seqs = degap.(seqs)
    resolve_seqs = resolve_seq.(seqs)
    kmer_vecs = kmer_count.(resolve_seqs,6)
    return collect(zip(seqnames,kmer_vecs,seqs))
end

function contam_check(filename,db; thresh = 0.015)
    seqnames,seqs = read_fasta_with_names(filename)
    seqs = degap.(seqs)
    resolve_seqs = resolve_seq.(seqs)
    kmer_vecs = kmer_count.(resolve_seqs,6)

    self_inds = [occursin(filename,s[1]) for s in db]
    if sum(self_inds) == 0
        @warn "No seqs from this sample in DB."
    end

    discard_calls = zeros(Bool,length(kmer_vecs))
    report_calls = zeros(Bool,length(kmer_vecs))
    closest_nonself_calls = ["" for i in 1:length(kmer_vecs)]
    closest_nonself_dists = zeros(length(kmer_vecs))

    for (i,kmer_vec) in enumerate(kmer_vecs)
        dists = [corrected_kmer_dist(kmer_vec,d[2]) for d in db]
        closest_self = minimum(dists[self_inds])
        closest_nonself = minimum(dists[.!self_inds])
        if closest_nonself < thresh
            report_calls[i] = true
            closest_nonself_calls[i] = (db[.!self_inds])[argmin(dists[.!self_inds])][1]
            closest_nonself_dists[i] = closest_nonself
            if closest_self > thresh
                discard_calls[i] = true
            end
        end
    end
    merged_calls = discard_calls .| report_calls
    #
    return seqnames[merged_calls],closest_nonself_calls[merged_calls],
      closest_nonself_dists[merged_calls],report_calls[merged_calls],discard_calls[merged_calls],discard_calls
end

