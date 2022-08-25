using DataFrames, NextGenSeqUtils, FASTX, WebBlast

"""
simple Hamming distance between two sequences
"""
function Hamming(s1, s2)
    if length(s1) != length(s2)
        @warn "cant compute Hamming dist on strings of unequal length"
        return max(length(s1),length(s2))
    end
    s1v=split(s1,"")
    s2v=split(s2,"")
    sgv=["-" for i in 1:length(s1)]
    return sum( (s1v .!= s2v) .& (s1v .!= sgv) .& (s2v .!= sgv) )
end

"""
compute pairwise distances between two sets of sequences
"""
function PairWise(f,xs,ys)
    return [f(x,y) for x in xs, y in ys]
end

"""
get Blast results for query sequences
"""
function get_blast_results(files; thresh_hold=0.05, max_clades=5, max_waits=30, verbose=false)
    queries = DataFrame(sample=String[], count=Int[], query=String[])
    for file in files
        sample = replace(basename(file), r".fasta" => "")
        ids, seqs = read_fasta(file)
        if endswith(file,"rejected.fasta")
            ms, msi = findmax(length.(degap.(seqs)))
            println("adding longest reject (sequence $(msi)), to blast queries...")
            ids = [ids[msi]]
            seqs = [seqs[msi]]
        end
        df = DataFrame(id=ids, seq=seqs)
        gdf = DataFramesMeta.groupby(df,:seq)
        gdf = combine(gdf,nrow => :count)
        sort!(gdf,:count,rev=true)
        cdf = similar(gdf,1)
        cdf[1,:] = gdf[1,:]
        gdf = gdf[2:end,:]
        clade = 1
        while nrow(gdf) > 0  && clade < max_clades
            dist_cols = PairWise(Hamming,gdf[!,:seq],cdf[!,:seq])
            dist_min = vec(mapslices(minimum,dist_cols,dims=2))
            edf = DataFrame(dm=dist_min)
            gdf = hcat(gdf,edf)
            sort!(gdf,:dm,rev=true)
            max_min_dist = gdf[1,:dm] / length(gdf[1,:seq])
            if verbose println("$(sample), max_min dist = ", max_min_dist) end
            if max_min_dist < thresh_hold
                if verbose println("$(sample), no more clades, breaking ...") end
                break
            end
            select!(gdf,Not(:dm))
            push!(cdf,gdf[1,:])
            gdf = gdf[2:end,[:seq,:count]]
            clade += 1
        end
        if verbose println("$(sample), $(clade) clade identified ...") end
        for i in 1:nrow(cdf)
            push!(queries, [sample*"_$(i)"
                            cdf[i,:count]
                            degap(cdf[i,:seq])])
        end
    end
    total_queries = nrow(queries)
    sample = replace(basename(files[1]), r".fasta" => "")
    println("$(sample), blasting $(total_queries) query sequences ... ")
    results = DataFrame(sample=String[],score=Float64[],accession=String[],description=String[])
    blast_arr = WebBLAST(queries[:,:query],query_names=queries[:,:sample],
                            max_waits=max_waits*total_queries, num_hits=1)
    if blast_arr != nothing
        rdf = flatten_to_dataframe(blast_arr)
        for (sample, identity, query_len, hit_acc, hit_def) in zip(rdf[:,"Query-def"], rdf[:,:Hsp_identity], rdf[:,"Query-len"], rdf[:,:Hit_accession], rdf[:,:Hit_def] )
            score = round(tryparse(Int64, identity) / tryparse(Int64, query_len),digits=3)
            if verbose
                println("sample = $(sample), score = $score, acc = $(hit_acc)")
                println("$(hit_def)")
            end
            push!(results, [sample score hit_acc hit_def])
        end
    end
    results = results[results[:,:score] .> 0.75, :]
    return(results)
end
