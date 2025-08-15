using Pkg
Pkg.activate("./")
Pkg.instantiate()
using BioSequences, FASTX, YAML, CSV, DataFrames
using StatsBase
using MAFFT_jll, FastTree_jll
using Plots, RecipesBase
using PhyloNetworks

function consensus(seqs)
    cons = join([mode([seqs[i][j]
                    for i in 1:length(seqs)])
                        for j in 1:length(seqs[1])])
    return(cons)
end

function mergeInput(in_file_1,in_file_2,out_file)
    d1, d2 = Dict(), Dict()
    stream = open(in_file_1)
    records = collect(FASTX.FASTA.Reader(stream))
    seqs = (ungap).( (x->FASTX.sequence(LongSequence{DNAAlphabet{4}},x)).(records) )
    nams = String.(FASTX.description.(records))
    umis = (x->split(x," ")[1][end-7:end]).(nams)
    for i in 1:length(umis)
        d1[umis[i]] = seqs[i]
    end
    weights = (x->parse(Int, split(x," ")[2][4:end])).(nams)
    seq_map_1=countmap(umis,weights)
    close(stream)
    stream = open(in_file_2)
    records = collect(FASTX.FASTA.Reader(stream))
    seqs = (ungap).( (x->FASTX.sequence(LongSequence{DNAAlphabet{4}},x)).(records) )
    nams = String.(FASTX.description.(records))
    umis = (x->split(x," ")[1][end-7:end]).(nams)
    for i in 1:length(umis)
        d2[umis[i]] = seqs[i]
    end
    weights = (x->parse(Int, split(x," ")[2][4:end])).(nams)
    seq_map_2=countmap(umis,weights)
    close(stream)
    both=intersect(keys(seq_map_1),keys(seq_map_2))
    @show length(both)
    only_1=setdiff(keys(seq_map_1),keys(seq_map_2))
    @show length(only_1)
    only_2=setdiff(keys(seq_map_2),keys(seq_map_1))
    @show length(only_2)
    stream = open(FASTA.Writer, out_file, append=false)
    count=0
    for k in both
        count+=1
        nam="seq$(count)_$(k)_$(seq_map_1[k])_$(seq_map_2[k])"
        write(stream, FASTA.Record(nam, d1[k]))
        if d1[k] != d2[k]
            @warn "mismatch at umi $(k)"
            count+=1
            nam="seq$(count)_$(k)_$(seq_map_1[k])_$(seq_map_2[k])"
            write(stream, FASTA.Record(nam, d2[k]))
        end
    end
    for k in only_1
        count+=1
        nam="seq$(count)_$(k)_$(seq_map_1[k])_0"
        write(stream, FASTA.Record(nam, d1[k]))
    end
    for k in only_2
        count+=1
        nam="seq$(count)_$(k)_0_$(seq_map_2[k])"
        write(stream, FASTA.Record(nam, d2[k]))
    end
    close(stream)
    return
end

function newickTree(seq_file, tree_file; nu=true)
    cmd=`$(fasttree_double()) -quiet -nosupport -out $(tree_file) $(seq_file)`
    if nu
        cmd=`$(fasttree_double()) -quiet -nt -gtr -nosupport -out $(tree_file) $(seq_file)`
    end
    cmd=`$fasttree -quiet -nt -gtr -nosupport -out $(tree_file) $(seq_file)`
    run(cmd);
    println("tree generated in $(tree_file)")
    return true
end

function my_mafft(in_file, out_file)
    arguments = "$(in_file) > $(out_file) "
    # run(`$(mafft_fftns()) -help`)
    # run(`$(mafft_fftns()) $(in_file) '>' $(out_file)`)
    # run(`mafft $(in_file) stdout = $(out_file)`)
    cmd = `mafft-fftns --quiet --thread 2 --ep 2 --op 3 --out $out_file $in_file`
    println(cmd)
    run(cmd)
    return
end

function find_most_common_seqname(seqnames)
    counts = (x->parse(Int, split(x,"_")[3])+parse(Int, split(x,"_")[4])).(seqnames)
    return seqnames[findmax(counts)[2]]
end

# reroot needs PhyloNetworks, try to do this with Phylo
function reroot(tree_file, rooted_file)
    tree = readTopology(read(tree_file,String))
    cons_name=find_most_common_seqname(tree.names) # "consensus"
    rootatnode!(tree, cons_name)
    directEdges!(tree)
    cladewiseorder!(tree)
    ladderize!(tree,getroot(tree))
    writeTopology(tree, rooted_file) 
    return(tree)
end
        
function get_max_depth(node)
    if isleaf(node)
        return getparentedge(node).length
    else
        return maximum( getparentedge(node).length .+ get_max_depth.(getchildren(node)) )
    end
end

function ladderize!(net,node)
    if isleaf(node)
        return 
    else
        children=getchildren(node)
        # sort the children on depths
        depths = get_max_depth.(children)
        if ! issorted(depths,rev=true)
            child_edge_numbers=(x->PhyloNetworks.getconnectingedge(node, x).number).(children)
            sp = sortperm(depths,rev=true)
            oens=child_edge_numbers[sp]
            PhyloNetworks.rotate!(net,node.number,orderedEdgeNum=oens)
        end
        (x->ladderize!(net,x)).(children)
        return
    end
end

function count_split(s,c)
    v=split(s,c)
    tc=(1.0+parse(Int,v[3])+parse(Int,v[4]))
    return tc
end

function getnodeheight(node,tree)
    if isrootof(node,tree)
        return 0.0
    else
        return getnodeheight(getparent(node),tree) + getparentedge(node).length
    end
end

function depth_first_heights(node,height,leaf_num)
    if isleaf(node)
        leaf_num += 1
        height[node.number]=leaf_num
        return leaf_num
    else
        children=getchildren(node)
        for child in children
            leaf_num=depth_first_heights(child,height,leaf_num)
        end
        height[node.number]=sum((x->height[x.number]).(children))/length(children)
        return leaf_num
    end
end

function get_color(name)
    parts = split(name,"_")
    if length(parts) == 4
        c1=parse(Float64,parts[3])
        c2=parse(Float64,parts[4])
        return c1/(c1+c2)
    else
        return nothing
    end
end
    
function findxy(tree::HybridNetwork)
    height, depth, names, colors = Dict(), Dict(), Dict(), Dict()
    heights=getnodeheights(tree,true)
    i=0
    for n in tree.node
        i=i+1
        height[n.number] = getnodeheight(n,tree)
        names[n.number] = n.name
        colors[n.number] = get_color(n.name)
    end
    leaf_num=0
    depth_first_heights(getroot(tree),depth,leaf_num)
    return height, depth, names, colors
end

@recipe function f(tree::HybridNetwork;
                    # nam_ali::Base.Iterators.Zip{Tuple{Vector{Any}, Vector{Any}}};
                    # ref_ali::Base.Iterators.Zip{Tuple{Vector{String}, Vector{String}}},
                    # ept::Dict{Any,Any};
                    label_top = "top",
                    label_bot = "bot",
                    treetype = :dendrogramhighlighter,
                    marker_group = nothing, line_group = nothing,
                    showtips = true, showtipmarkers = true, tipfont = (4,"Courier Bold"))
    
    linecolor --> :black
    linewidth --> 3
    grid --> false
    framestyle --> :none
    legend --> false
    colorbar --> true
    # size --> (1000, 1000)

    tip_numbers = (x->x.number).(tree.leaf)
    tip_names = (x->x.name).(tree.leaf)
    
    d, h, n, c = findxy(tree)
    
    size_x=80*14
    size_y=(length(tip_numbers))*10
    @show maximum(values(d))
    scale_x=size_x / (maximum(values(d)))
    @show size_x, size_y, scale_x
    for k in keys(d)
        d[k]=scale_x*d[k]
    end
    @show maximum(values(d))
    size --> ( size_x , size_y  )
    
    adj = 2 * 14 # 0.1 * maximum(values(d))
    @show adj
    adj <= 0.0 ? adj=0.1 : nothing
    
    tipmarker_x = map(x -> d[x], tip_numbers)
    tipmarker_y = map(x -> h[x], tip_numbers)
    tipmarker_s = map(x -> log10(10+count_split(x,'_')), tiplabels(tree))
    tipmarker_g = map(x -> c[x], tip_numbers)
    # global visits = vcat(["consensus"],sort(union(tipmarker_g))[1:end-1])
    max_x = maximum(values(d))
    min_y = minimum(values(h))
    # anot_dic, ept_anot = getNameAnnotation(ept, ref_ali, nam_ali)
    # tipannotations = map(x -> (max_x + adj, h[x], anot_dic[n[x]]), tip_numbers)
    tipannotations = map(x -> (d[x] + 14, h[x], n[x], c[x]), tip_numbers)  # * log10(3+count_split(n[x],'_'))
    

    x, y = Float64[], Float64[]
    
    for node in tree.node
        try p = getparent(node) 
            push!(x, d[p.number], d[p.number],
                  d[node.number], NaN)
            push!(y, h[p.number], h[node.number], h[node.number], NaN)
        catch e
        end
    end
    # add the scale bar
    bar=0.01*scale_x
    push!(x,0,bar,NaN,0,0,NaN,bar,bar,NaN)
    push!(y,0,0,NaN,-0.5,0.5,NaN,-0.5,0.5,NaN)
    push!(tipannotations,(bar+14,0,"0.01",0.0))

    marker_x, marker_y, marker_z = values(d), values(h), values(c)

    
    if treetype == :dendrogramhighlighter
        DendrogramHighlighter(x, y, label_top, label_bot, tipannotations, marker_x, marker_y, marker_z,
                    tipmarker_x, tipmarker_y, tipmarker_s, tipmarker_g,
                    showtips, showtipmarkers, tipfont)
    elseif treetype == :fanhighlighter
        FanHighlighter(x, y, label_top, label_bot, tipannotations, marker_x, marker_y, marker_z,
                    tipmarker_x, tipmarker_y, tipmarker_s, tipmarker_g,
                    showtips, showtipmarkers, tipfont)
    else
        throw(ArgumentError("Unsupported `treetype`; valid values are `:dendrogramhighlighter` or `:fanhighlighter`"))
    end
end

struct DendrogramHighlighter
    x::Any
    y::Any
    label_top::Any
    label_bot::Any
    tipannotations::Any
    marker_x::Any
    marker_y::Any
    marker_z::Any
    tipmarker_x::Any
    tipmarker_y::Any
    tipmarker_s::Any
    tipmarker_g::Any
    showtips::Any
    showtipmarkers::Any
    tipfont::Any
    # marker_group::Any
    # line_group::Any
end

struct FanHighlighter
    x::Any
    y::Any
    label_top::Any
    label_bot::Any
    tipannotations::Any
    tipmarkers::Any
    marker_x::Any
    marker_y::Any
    marker_z::Any
    tipmarker_x::Any
    tipmarker_y::Any
    tipmarker_s::Any
    tipmarker_g::Any
    showtips::Any
    showtipmarkers::Any
    tipfont::Any
    # marker_group::Any
    # line_group::Any
end

@recipe function f(dend::DendrogramHighlighter)
    ex = extrema(filter(isfinite, dend.x))
    # xlims --> (ex[1] - 0.5 * ex[2], ex[2] * 6.0)
    xlims --> (ex[1] - 5*14 , ex[2] + 25*14)
        
    ey = extrema(filter(isfinite, dend.y))
    ylims --> (ey[1] - 2.0, ey[2] + 4.0)
    
   

    # tip annotations
    @show dend.tipfont
    
    grads = cgrad(:darkrainbow)

    dend.showtips &&
        ( annotations := vcat( map(x -> (x[1], x[2],
                            text(x[3], grads[Int(max(1,ceil(x[4]*length(grads))))], :left, dend.tipfont...)),
                            dend.tipannotations) ,
                            [ ( ex[2] + 30*8, ey[2]+1, text(dend.label_top, grads[end], :left, 8, "Courier Bold") ),
                              ( ex[2] + 30*8, ey[1]-1, text(dend.label_bot, grads[1], :left, 8, "Courier Bold") ) ]
                            ) )
    
    @series begin
        seriestype := :path
        markersize := 0
        markershape := :none
        series_annotations := nothing
        label --> "path"
        dend.x, dend.y
    end

    if dend.showtipmarkers
        @series begin
            seriestype := :scatter
            label --> "mix"
            clims := (0,1)
            c := :darkrainbow
            colorbar_ticks := (0:1.0:1.0, string.(round.(Int, (0:1.0:1.0) .* 100), "%"))
            markersize := (x->x).(dend.tipmarker_s)
            marker_z := (x->x).(dend.tipmarker_g)
            zcolor := range(0.0,stop=1.0,length=5)
            markerstrokewidth := 0
            # markercolor := (x->weighted_color_mean(x, colorant"red", colorant"green")).(dend.tipmarker_g)
            # markerstrokecolor := (x->weighted_color_mean(x, colorant"red", colorant"green")).(dend.tipmarker_g)
            dend.tipmarker_x, dend.tipmarker_y
        end
    end
            
    primary := false
    label := "something"
    return nothing
end

####################### script starts now ###################################

if length(ARGS) != 2 || endswith(".fasta", ARGS[1]) && endswith(".fasta", ARGS[2])
    println("*******************************************************************")
    println("usage: julia MergeTreeplot.jl pacbio_fasta_file nanopore_fasta_file")
    println("*******************************************************************")
    exit()
end

in_file_1=ARGS[1]
in_file_2=ARGS[2]

work_dir = "working/"
mkpath(work_dir)
merged_file = work_dir*"merged.fasta"
merged_aligned_file = work_dir*"merged_aligned.fasta"
tree_file = work_dir*"tree.tre"
rerooted_tree_file = work_dir*"rerooted_tree.tre"
ladder_tree_file = work_dir*"ladder_tree.tre"
plot_file = work_dir*"fs_tree_plot_$(basename(in_file_1)[1:end-6])_$(basename(in_file_2)[1:end-6]).pdf"

mergeInput(in_file_1,in_file_2,merged_file)
my_mafft(merged_file, merged_aligned_file)
newickTree(merged_aligned_file, tree_file; nu=false)

tree=reroot(tree_file,rerooted_tree_file)
PhyloNetworks.resetnodenumbers!(tree; checkpreorder=true, type=:postorder)
ladderize!(tree,getroot(tree))

writeTopology(tree, ladder_tree_file)

title="\n Merged UMIs \n $(basename(in_file_1)[1:end-6]) \n $(basename(in_file_2)[1:end-6])"

stream = open(merged_aligned_file)
records = collect(FASTX.FASTA.Reader(stream))
seqs = (ungap).( (x->FASTX.sequence(LongSequence{DNAAlphabet{4}},x)).(records) )
nams = String.(FASTX.description.(records))
close(stream)

nam_ali = zip(nams,seqs)

pl=plot(tree, linecolor = :orange, linewidth = 1,
        label_top = split(basename(in_file_1),"_")[1],
        label_bot = split(basename(in_file_2),"_")[1],
        showtips = true, showtipmarkers=true, aligntips=true,
        treetype = :dendrogramhighlighter,
        title=title, titlefontsize=8 )
        
savefig(pl,plot_file)
println("Figure saved in $(plot_file)" )
exit()
