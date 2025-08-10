using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using CSV, DataFrames, Plots
using BioSequences, FASTX, YAML, StatsBase
using CodecZlib: GzipDecompressorStream
using CodecZlib: GzipCompressorStream

################################## some functions #########################

function string_reverse_complement(dna_string::String)
    return String(reverse_complement(LongDNA{4}(dna_string)))
end

function primer_peek(seqs::Vector{Any}; L=20, N=30, keep=20 )
    starts=[ s[1:L] for s in seqs ]
    ends=[ s[end-L+1:end] for s in seqs ]
    pot_primers = countmap( vcat(starts,string_reverse_complement.(ends)))
    topN=reverse(sort([(pot_primers[k],k) for k in keys(pot_primers)]))[1:N]
    for i in enumerate(topN)
        println(i)
    end
    # return my_reverse_complement.(sort(my_reverse_complement.([i[2] for i in topN[1:keep]])))
    return topN[1:keep]
end

function findbestlast(query,tol,window)
    best=nothing
    next=findfirst(query,tol,window)
    while !isnothing(next)
        best=next
        next=findnext(query,tol,window,last(next)+1)
    end
    return best
end



##################### main script starts here ##############################

if length(ARGS) < 2
    println("******************************************************")
    println("usage: julia umi_analysis.jl pacbio_file nanopore_file")
    println("******************************************************")
    exit()
end

pacbio_file=ARGS[1]
nanopore_file=ARGS[2]

t1 = time()

pb_df=CSV.read(pacbio_file,DataFrame)
pb_df=pb_df[pb_df[!,:tags].=="fs<5" .|| pb_df[!,:tags].=="likely_real",:]
np_df=CSV.read(nanopore_file,DataFrame)
np_df=np_df[np_df[!,:tags].=="fs<5" .|| np_df[!,:tags].=="likely_real",:]
@show size(pb_df)
@show size(np_df)

pb_umis=setdiff(pb_df[!,:UMI],np_df[!,:UMI])
np_umis=setdiff(np_df[!,:UMI],pb_df[!,:UMI])
pb_only=pb_df[(x->in(x,pb_umis)).(pb_df[!,:UMI]).&&(pb_df[!,:fs].>=5),:]
np_only=np_df[(x->in(x,np_umis)).(np_df[!,:UMI]).&&(np_df[!,:fs].>=5),:]

# @show pb_umis
# @show pb_only[:,:fs]
# @show np_umis
# @show np_only[:,:fs]
# histogram(pb_only[:,:fs],label="PacBio",alpha=0.4)
# histogram!(np_only[:,:fs],label="Nanopore",alpha=0.4)
# savefig("$(split(basename(pacbio_file),'_')[2])_fs_not_common.png")

pb_jitter_x=rand(size(pb_only)[1]).*0
pb_jitter_y=rand(size(pb_only)[1]).*0
np_jitter_x=rand(size(np_only)[1]).*0
np_jitter_y=rand(size(np_only)[1]).*0

common_umis=intersect(pb_df[!,:UMI],np_df[!,:UMI])


pb_df=pb_df[(x->in(x,common_umis)).(pb_df[!,:UMI]),:]
np_df=np_df[(x->in(x,common_umis)).(np_df[!,:UMI]),:]
@show size(pb_df)
@show size(np_df)

# rename!(pb_df,:fs => :pb_fs)
# rename!(np_df,:fs => :np_fs)

df = innerjoin(pb_df,np_df,on = :UMI, makeunique=true)

scatter(df[:,:fs],df[:,:fs_1],xlabel="PacBio fs",ylabel="Nanopore fs",label="PacBio and Nanopore",alpha=0.3,
            title="$(split(basename(pacbio_file),'_')[2]): Family Sizes for Commom UMIs")
scatter!(pb_only[:,:fs].+pb_jitter_x,pb_only[:,:fs].+pb_jitter_y,label="PacBio only",alpha=0.9,c="red")
scatter!(np_only[:,:fs].+np_jitter_x,np_only[:,:fs].+np_jitter_y,label="Nanopore only",alpha=0.9,c="orange")
savefig("$(split(basename(pacbio_file),'_')[2])_fs_common.png")

# now do a minag scatterplot
pb_df=CSV.read(pacbio_file,DataFrame)
pb_df=pb_df[pb_df[!,:tags].=="fs<5" .|| pb_df[!,:minag].>0.0,:]
np_df=CSV.read(nanopore_file,DataFrame)
np_df=np_df[np_df[!,:tags].=="fs<5" .|| np_df[!,:minag].>0.0,:]

common_umis=intersect(pb_df[!,:UMI],np_df[!,:UMI])
pb_df=pb_df[(x->in(x,common_umis)).(pb_df[!,:UMI]),:]
np_df=np_df[(x->in(x,common_umis)).(np_df[!,:UMI]),:]
df = innerjoin(pb_df,np_df,on = :UMI, makeunique=true)

df=df[ (df[:,:minag].>0.0) .&& (df[:,:minag_1].>0.0),:]
scatter(df[:,:minag],df[:,:minag_1],xlabel="PacBio minag",ylabel="Nanopore minag",
    label="common PacBio and Nanopore UMIs",alpha=0.3,xlim=[0,1],ylim=[0,1],
    title="$(split(basename(pacbio_file),'_')[2]): Minimum Agreement for Commom UMIs")
vline!([0.7],c="orange",label="Existing PacBio threshold 0.7")
hline!([0.6],c="red",label="Recommended Nanopore threshold 0.6")
savefig("$(split(basename(pacbio_file),'_')[2])_minag_common.png")

t2 = time()

dt=Int( floor( (t2-t1)/60 ) )

exit()
