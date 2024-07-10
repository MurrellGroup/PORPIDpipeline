using NextGenSeqUtils, DPMeansClustering,
RobustAmpliconDenoising, StatsBase,
PyPlot, MultivariateStats,
BioSequences,
Distributions,LinearAlgebra,DataFrames,CSV,
DataFramesMeta, Seaborn,
Compose, Colors

import MolecularEvolution

#Extra import
function read_fasta_with_everything(filename; seqtype=String)
    records = NextGenSeqUtils.read_fasta_records(filename)
    return seqtype[FASTA.sequence(seqtype, r) for r in records], [FASTA.identifier(r) for r in records], [FASTA.description(r) for r in records]
end

function parse_HVTN_descriptor(descriptor)
    try
        splot = split(descriptor,"=")
        return parse(Int64,split(splot[2]," ")[1]),parse(Float64,splot[3])
    catch
        @warn "Broken descriptor"
        return -999,-999.0
    end
end

function read_fasta_with_descriptors_in_names(f)
   seqs,seqnames,descriptors = read_fasta_with_everything(f);
   seqnames = seqnames .* " " .* descriptors
   return seqnames,seqs
end

#Various distances.
function non_gap_hamming(s1,s2)
    l = min(length(s1),length(s2))
    diff = 0
    for i in 1:l
        if s1[i] != '-' && s2[i] != '-'
            if s1[i] != s2[i]
                diff += 1
            end
        end
    end
    return diff
end

function normalized_non_gap_hamming(s1,s2)
   l = min(length(s1),length(s2))
   diff = 0
   count = 0
   for i in 1:l
       if s1[i] != '-' && s2[i] != '-'
           count += 1
           if s1[i] != s2[i]
               diff += 1
           end
       end
   end
   return diff/count
end

function pairwise_dist(seqs,dist_func)
    l = length(seqs)
    mat = zeros(l,l)
    for i in 1:l
        for j in i+1:l
            d = dist_func(seqs[i],seqs[j])
            #Assumes the distance is symmetric
            mat[i,j] = d
            mat[j,i] = d
        end
    end
    return mat
end

"""
    mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
Julia wrapper for mafft.
"""
function mafft(inpath, outpath; path="", flags::Vector{String}=String[], kwargs...)
    mktempdir() do mydir
       mafft = length(path) == 0 ? PATHS.mafft : path
       flagstrings = []
       for flag in flags
          push!(flagstrings, "--$flag")
       end
       args = []
       for (k, v) in kwargs
          push!(args, "--$k")
          push!(args, "$v")
       end
       progress = "$(mydir)/mafft.progress"
       # --adjustdirection
       run(`$mafft $flagstrings $args --quiet --progress $progress --out $outpath $inpath`)
   end
end

"""
    mafft_align{T<:BioSequence}(seqs::Vector{T}; kwargs...)
Julia wrapper for mafft.
"""
function mafft_align(seqs; kwargs...)
    mktempdir() do mydir
        seqfile = string(mydir, "/sequences.fasta")
        mafftout = string(mydir, "/mafft.fasta")
        write_fasta(string(mydir, "/sequences.fasta"), seqs)
        mafft(seqfile, mafftout; kwargs...)
        aligned = Array{String}(read_fasta(mafftout)[2])
        if  sum(uppercase.(degap.(seqs)) .!= uppercase.(degap.(aligned))) > 0
            println(" error in mafft, degapped input not equal to degapped output")
            in_size = length(seqs)
            out_size = length(aligned)
            if in_size == out_size
                for k in 1:in_size
                    if  uppercase.(degap.(seqs[k])) != uppercase.(degap.(aligned[k]))
                        println("seq_no = $(k)")
                        println(uppercase.(degap.(seqs[k])))
                        println(uppercase.(degap.(aligned[k])))
                    end
                end
            end
            @error "Aligned seqs don't match input seqs (input has $(in_size) seqs, output has $(out_size) seqs)!"
            exit(1)
        end
        return aligned
    end
end

function gap_preserving_consensus(ali_seqs)
    charvec = collect.(ali_seqs)
    return join([mode([c[i] for c in charvec]) for i in 1:length(charvec[1])])
end

function pairwise_G2A(s1,s2)
    l = min(length(s1),length(s2))
    diff = 0
    GA = 0
    for i in 1:l
        if s1[i] != '-' && s2[i] != '-'
            if s1[i] != s2[i]
                diff += 1
            end
            if (s1[i] == 'g' || s1[i] == 'G') && (s2[i] == 'a' || s2[i] == 'A')
                GA += 1
            end
        end
    end
    if diff < 1
        return 0.0
    else
        return GA/(diff + 20)
    end
end

function G2A_ratio(ali_seqs)
    cons = gap_preserving_consensus(ali_seqs)
    return pairwise_G2A.(cons,ali_seqs)
end

function MDS_plot(mds;
        sizes = ones(length(seqs)),
        color_vec = ones(length(seqs)),
        figsize = (6,5),
        default_sizes = [10,50,100],
        color_limits = (0.5,1.0),
        cmap = "jet",
        alpha = 0.7
    )
    fig = figure(figsize=figsize)
    for d in default_sizes
        scatter([0.0],[-99999.0],s = [d],color = "black",label = string(d),alpha = alpha)
    end
    init_xspan = maximum(mds[1,:]) - minimum(mds[1,:])
    init_yspan = (maximum(mds[2,:]) - minimum(mds[2,:]))

    if init_yspan > init_xspan
        mds[1,:],mds[2,:] = mds[2,:],mds[1,:]
    end
    scatter(mds[1,:],mds[2,:],s=sizes,cmap = cmap,c = color_vec, alpha = alpha)
    clim(color_limits[1], color_limits[2])
    starty = (minimum(mds[2,:]),maximum(mds[2,:]))
    x_span = maximum(mds[1,:]) - minimum(mds[1,:])
    padding = (x_span - (starty[2]-starty[1]))/2
    xlim((minimum(mds[1,:])-(x_span/10),maximum(mds[1,:])+(x_span/10)))
    ylim((starty[1]-padding-(x_span/10),starty[2]+padding+(x_span/10)))
    colorbar()
    #xlabel("MDS 1")
    #ylabel("MDS 2");
    legend()

    #println("xlims: $((minimum(mds[1,:]), maximum(mds[1,:])))")
    #println("ylims: $((starty[1]-padding,starty[2]+padding))")
    #println("xspan: $(x_span)")
    #println("yspan: $((starty[2]+padding)-(starty[1]-padding))")
    return fig
end

# Possible to generalize parsing of the descriptors?
function read_HVTN(f)
    seqs,seqnames,descriptors = read_fasta_with_everything(f);
    numeric_descriptors = parse_HVTN_descriptor.(descriptors)
    sizes = [i[1] for i in numeric_descriptors]
    agreement_scores = [i[2] for i in numeric_descriptors];
    return seqs,seqnames,sizes,agreement_scores
end

function align_HVTN(f)
    seqs,seqnames,sizes,agreement_scores = read_HVTN(f)
    ali_seqs = mafft_align(seqs);
    return ali_seqs,seqnames,sizes,agreement_scores
end

function fasta2aligned(f)
    @warn "This function the same as `align_HVTN(f)`"
    seqs,seqnames,sizes,agreement_scores = read_HVTN(f)
    ali_seqs = mafft_align(seqs);
    return ali_seqs,seqnames,sizes,agreement_scores
end

function ali2MDS(ali_seqs,sizes,colors) #peraps another MDS allowing annotations
    dist_mat = pairwise_dist(ali_seqs,normalized_non_gap_hamming);
    mds = MultivariateStats.transform(fit(MDS, dist_mat, maxoutdim=2, distances=true))
    return MDS_plot(mds,sizes = sizes, color_vec = colors)
end

function extract_region(ali_seqs,ref_profile)
    sample_profile = seqs2profile(uppercase.(ali_seqs));
    ali_ref, ali_sample = profile_affine_align(ref_profile, sample_profile,profile_cost);
    matched_inds = [i[1][1] for i in ali_sample] .!= '!';
    region_inds = [i[1][1] for i in ali_ref] .!= '!';
    region_start,region_end = findfirst(region_inds[matched_inds]),findlast(region_inds[matched_inds]);
    extracted_seqs = [s[region_start:region_end] for s in ali_seqs]
    return extracted_seqs,(region_start,region_end)
end

function extract_and_score_misalignments(ali_seqs,ref_profile; prob_func = get_prob_ignoring_indels)
    sample_profile = seqs2profile(uppercase.(ali_seqs));
    ali_ref, ali_sample = profile_affine_align(ref_profile, sample_profile,profile_cost);
    matched_inds = [i[1][1] for i in ali_sample] .!= '!';
    region_inds = [i[1][1] for i in ali_ref] .!= '!';
    region_start,region_end = findfirst(region_inds[matched_inds]),findlast(region_inds[matched_inds]);
    extracted_seqs = [s[region_start:region_end] for s in ali_seqs]

    region_probs = ali_ref[matched_inds][region_start:region_end]
    scores = [i[1] for i in score_sequence.(uppercase.(extracted_seqs);
            matched_ref_probs = region_probs,
            prob_func = get_prob_ignoring_indels)]
    return extracted_seqs, scores
end

function get_prob(char, prob_dict)
    if char in keys(prob_dict) #returns prob for observed
        prob = prob_dict[char]
    else
        prob = 0. #returns zero if not observed in panel
    end
    return prob
end

function get_prob_ignoring_indels(char, prob_dict)
    if char == '-' #ignores deletions
        prob = nothing
    elseif '!' in keys(prob_dict) #ignores insertions relative to panel
        prob = nothing
    elseif char in keys(prob_dict) #returns prob for observed
        prob = prob_dict[char]
    else
        prob = 0. #returns zero if not observed in panel
    end
    return prob
end

function get_prob_ignoring_ins(char, prob_dict)
    if '!' in keys(prob_dict) #ignores insertions relative to panel
        prob = nothing
    elseif char in keys(prob_dict) #returns prob for observed
        prob = prob_dict[char]
    else
        prob = 0. #returns zero if not observed in panel
    end
    return prob
end

#Kadane's algorithms
function max_subarray(arr::Vector{<:Number})
   largest_ending_here = 0
   best_start = this_start = ending = best_so_far = 1
   for (i,x) in enumerate(arr)
       largest_ending_here += x
       best_so_far = max( best_so_far, largest_ending_here )
       if largest_ending_here <= 0
           this_start = i + 1
           largest_ending_here = 0
       elseif largest_ending_here == best_so_far
           best_start = this_start
           ending = i
       end
   end
   return best_so_far, best_start, ending
end

function score_sequence(seq;
    prob_func = get_prob_ignoring_indels, matched_ref_probs = matched_ref_probs)
    probs = prob_func.([c for c in seq],[Dict(i) for i in matched_ref_probs]);
    probs = [p for p in probs if !isnothing(p)];
    probs_trans = -log.(probs .+ 10^-2) .+ log(1/4); #tweak scoring here
    return max_subarray(probs_trans)
end

function robust_translate(s)
    s = replace(s,"-"=>"")
    s = s[1:3*Int64(floor(length(s)/3))]
    translate_to_aa(s)
end

function find_first_stop(AAseq)
    s = collect(AAseq)
    if '*' in s
        return findfirst(s .== '*')
    else
        return length(AAseq)
    end
end

function next_frame_ct_gaps(s;ix=0)
    ct = 0
    frame =""
    for i in 1:length(s)
        ix+=1
        if s[i] != '-'
            ct+=1
            frame=frame*"$(s[i])"
            if ct == 3
            return frame,ix,s[i+1:end]
            end
        end
    end
    return "",ix,""
end

function find_first_stop_pos_with_gaps(s;ix=0)
    frame, ix, s_j = next_frame_ct_gaps(s;ix=ix)
    #println(frame)
    if frame in ["TAA","TAG","","TGA"]
        return ix
    else
        find_first_stop_pos_with_gaps(s_j;ix=ix)
    end
end

function env_function_check_2(AAseq; min_ORF_length = 820)
    first_stop = find_first_stop(AAseq)
    flag = true
    reason = ""
    if first_stop < min_ORF_length
        reason *= "Short ($(length(AAseq)) ORF) - deletion or early stop."
        flag = false
    end
    return flag, AAseq[1:first_stop], reason
end

#Returns: flag, seq, reason
function env_function_check(AAseq; min_degapped_AA_length = 820, max_stop_dist_from_end = 150)
    #Note: This function will discard sequences that are too short because of a large deletion,
    #but keep those that end a little early due to a stop-induced CT truncation. If a deletion
    #occurs in the CT, these also get discarded, which is not ideal.
    first_stop = find_first_stop(AAseq)
    flag = true
    reason = ""
    if length(AAseq) < min_degapped_AA_length
        reason *= "Short ($(length(AAseq)) AAs including stops) - likely deletion. "
        flag = false
    end
    if length(AAseq) - first_stop > max_stop_dist_from_end
        reason *= "Stop $(length(AAseq) - first_stop) from EOS. "
        flag = false
    end
    return flag, AAseq[1:first_stop], reason
end

function mafft_stop_preserving_align(AA_seqs)
    modified_seqs = [replace(s,"*"=>"B") for s in AA_seqs]
    AA_ali = mafft_align(modified_seqs)
    replaced_seqs = [replace(s,"B"=>"*") for s in AA_ali]
    return replaced_seqs
end


function get_mut_counts(count_mat)
    n_mutations = sum(count_mat) - sum(diag(count_mat))
    return (count_mat[3,1], n_mutations)
end

#=
function ali2MDS_anno(ali_seqs,sizes,agreement_scores,labels; export_path = "out.pdf")
    dist_mat = pairwise_dist(ali_seqs,normalized_non_gap_hamming);
    mds = transform(fit(MDS, dist_mat, maxoutdim=2, distances=true))
    MDS_plot(mds,sizes = sizes, color_vec = agreement_scores)
    for (ix,label) in enumerate(labels)
        PyPlot.text(x = mds[1,ix], y = mds[2,ix]; s = label, fontsize = 8)
    end
    if !isnothing(export_path)
        savefig(export_path)
    end
end

=#

function variant_collapse(seqs; prefix = "seq_")
    dic = countmap(uppercase.(seqs))
    OC = sort(dic, by = x -> dic[x], rev = true)
    seqs,sizes = collect(keys(OC)),collect(values(OC))
    names = [prefix*string(i)*"_"*string(sizes[i]) for i in 1:length(seqs)]
    return seqs,sizes,names
end

#function collapse_nucs_by_AA(ali_nuc,ali_AA)
#    AA
#end

function fasttree_nuc(seqs,seqnames; quiet = false)
   mktempdir() do mydir
       seqfile = string(mydir, "/sequences.fasta")
       treefile = string(mydir, "/fasttree.newick")
       write_fasta(string(mydir, "/sequences.fasta"), seqs,names = seqnames)
       if quiet
       run(`fasttree -quiet -nosupport -nt -gtr -out $(treefile) $(seqfile)`)
       else
       run(`fasttree -nosupport -nt -gtr -out $(treefile) $(seqfile)`)
       end
       #open the file and return lines
       lines = open(string(mydir,"/fasttree.newick"), "r") do io
           readlines(io)
       end
       return lines
   end
end

function fasttree_AA(seqs,seqnames; quiet = false)
   mktempdir() do mydir
       seqfile = string(mydir, "/sequences.fasta")
       treefile = string(mydir, "/fasttree.newick")
       write_fasta(string(mydir, "/sequences.fasta"), seqs,names = seqnames)
       if quiet
       run(`fasttree -quiet -nosupport -out $(treefile) $(seqfile)`)
       else
       run(`fasttree -nosupport -out $(treefile) $(seqfile)`)
       end
       #open the file and return lines
       lines = open(string(mydir,"/fasttree.newick"), "r") do io
           readlines(io)
       end
       return lines
   end
end

"""
    artefact_cutoff(ccs_counts,af_thresh)
Compute ccs cutoff from vector of ccs_counts and a threshold between 0 and 1
"""
function artefact_cutoff(ccs_counts,af_thresh)
    # ccs=sort(ccs_counts)
    # tot=sum(ccs)
    # cum_ccs=cumsum(ccs)
    # cut=tot*af_thresh
    # cut_ind=findfirst(x->x>cut,cum_ccs)
    # af_cutoff=ccs[cut_ind]
    af_cutoff = Int(ceil(maximum(ccs_counts)*af_thresh))
    return(af_cutoff)
end

"""
    family_size_umi_len_stripplot
Draws a stripplot of family sizes vs. UMI length from an input
DataFrame. Returns the figure object.
"""
function family_size_umi_len_stripplot(data;
                    fs_thresh=5, af_thresh=0.15, af_cutoff=1)
    tight_layout()
    fig = figure(figsize = (6,2))
    ax = PyPlot.axes()

    stripplot(y = [length(ix) for ix in data[!,:UMI]],
        x = data[!,:fs],
        hue = data[!,:tags],
        hue_order = [
            "fs<$(fs_thresh)",
            "possible_artefact",
            "likely_real",
            "UMI_len != 8",
            "LDA-rejects",
            "heteroduplex"
        ],
        alpha = 0.2, dodge = true, jitter = 0.3, orient = "h")
    
    # af_cutoff=artefact_cutoff(data[!,:fs], af_thresh)
    axvline([af_cutoff-0.5],c="red",label="artefact threshold")
    
    aftp=Int(af_thresh*100)
    labels = xlabel("UMI family size"), ylabel("UMI length")
    #  with $(aftp)% artefact threshold (fs=$(af_cutoff))
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # Summary for plot title
    t = @transform(data, :is_likely_real = :tags .== "likely_real")
    g = DataFramesMeta.groupby(t,:is_likely_real)
    # cts = @based_on(g, CCS = sum(:fs))
    cts = @combine(g, :CCS = sum(:fs))
    cts = sort!(cts, [:is_likely_real], rev = true)
    title("likely_real: $(cts[1, :CCS]) CCS, rejected: $(cts[2, :CCS]) CCS")

    return fig
end

"""
    family_size_stripplot
Draws a stripplot of family sizes with large jitter for possible_artefacts
and likely_reals from an input DataFrame. Returns the figure object.
"""
function family_size_stripplot(data;
                    fs_thresh=5, af_thresh=0.15, af_cutoff=1)
    tight_layout()
    fig = figure(figsize = (6,2))
    ax = PyPlot.axes()

    stripplot( y = [length(ix) for ix in data[!,:UMI]],
        x = data[!,:fs],
        hue = data[!,:tags],
        hue_order = [
            "fs<$(fs_thresh)",
            "possible_artefact",
            "likely_real"
        ],
        alpha = 0.6, dodge = false, jitter = 0.8, orient = "h")
    
    # af_cutoff=artefact_cutoff(data[!,:fs], af_thresh)
    
    axvline([af_cutoff-0.5],c="red",label="artefact threshold")
    
    for afths in 0.05:0.1:0.85
        afc=artefact_cutoff(data[!,:fs], afths)
        axvline([afc-0.5],alpha=1.0)
    end
    
    afths=0.95
    afc=artefact_cutoff(data[!,:fs], afths)
    axvline([afc-0.5],alpha=1.0,label="5% 15% ... 95%")
    
    axvline([af_cutoff-0.5],c="red")
    
    
    aftp=Int(af_thresh*100)
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    # Summary for plot title
    # t = @transform(data, :is_likely_real = :tags .== "likely_real")
    # g = DataFramesMeta.groupby(t,:is_likely_real)
    # cts = @based_on(g, CCS = sum(:fs))
    # cts = @combine(g, :CCS = sum(:fs))
    # cts = sort!(cts, [:is_likely_real], rev = true)
    # pc_artefact = round(100 * cts[2, :CCS] / (cts[1, :CCS] + cts[2, :CCS]),digits=1)
    nlr=sum(data[!,:tags].=="likely_real")
    naf=sum(data[!,:tags].=="possible_artefact")
    pc_artefact = round(100 * naf / (nlr + naf),digits=1)
    
    labels = xlabel("$(aftp)% af-thresh (fs=$(af_cutoff)) = $(pc_artefact)% artefacts"),
            ylabel("jitter plot")
    # title("percent artefact = $(pc_artefact)%")

    return fig
end


"""
    di_nuc_freqs(umis,weights)

    Draws two annotated heatmaps of di-neucleotide frequencies in a list of strings,
    the first unweighted and the second weighted.
"""
function di_nuc_freqs(umis::Vector{<:String}; weights=ones(length(umis)))
    tight_layout()
    # ioff()
    fig = figure("dinucleotide_table",figsize = (8,3))
    # fig.suptitle(title)
    alphabet = join(union(join(umis)))
    n = length(alphabet)
    freq_mat = reshape(zeros(n*n),(n,n))
    freq_mat_weighted = reshape(zeros(n*n),(n,n))
    count = 0
    count_weights = 0
    for m in 1:length(umis)
      for k in 1:length(umis[m])-1
        i = findall(string(umis[m][k]), alphabet)[1][1]
        j = findall(string(umis[m][k+1]), alphabet)[1][1]
        freq_mat[i,j] += 1
        freq_mat_weighted[i,j] += weights[m]
        count += 1
        count_weights += weights[m]
      end
    end
    if count > 0
      freq_mat = freq_mat / count
    end
    freq_mat = reshape(vcat([ reverse(freq_mat[:,j]) for j in 1:n ]...),(n,n))
    subplot(121)
    p = heatmap(freq_mat,annot=true,square=true, fmt=".2f",
      vmin=0, vmax=0.2, cbar=true, cmap="YlOrBr",
      xticklabels=string.(union(alphabet)),
      yticklabels=string.(reverse(union(alphabet))),
      cbar_kws=Dict("ticks" => [0.0, 0.05, 0.1, 0.15, 0.2],
                    "format" => "%.2f")) # "label" => "UMI families"))
    p.set_title("UMI families")
    if count_weights > 0
      freq_mat_weighted = freq_mat_weighted / count_weights
    end
    freq_mat_weighted = reshape(vcat([ reverse(freq_mat_weighted[:,j]) for j in 1:n ]...),(n,n))
    subplot(122)
    p = heatmap(freq_mat_weighted,annot=true,square=true, fmt=".2f",
      vmin=0, vmax=0.2, cbar=true, cmap="YlOrBr",
      xticklabels=string.(union(alphabet)),
      yticklabels=string.(reverse(union(alphabet))),
      cbar_kws=Dict("ticks" => [0.0, 0.05, 0.1, 0.15, 0.2],
                    "format" => "%.2f")) #"label" => "all UMIs"))
    p.set_title("all UMIs")
    close(fig)
    return fig
end

function leaf_rename(str)
    #s1 = split(str,"H704_")
    return str #s1[1]*join(split(s1[2],"_")[5:6],"_")
end

function get_size(s)
    return parse(Int64,split(s,"_")[end])
end

function draw_variant_tree(treestring, rename_func, filepath)
    newt = gettreefromnewick(treestring, GeneralFelNode);
    max_dist = maximum(root2tip_distances(newt)[1])
    newt.name = string(round(max_dist,sigdigits = 4))

    for n in getleaflist(newt)
        n.name = rename_func(n.name)
    end

    dot_size_dict = Dict()
    all_sizes = []
    for n in getleaflist(newt)
        dot_size_dict[n] = sqrt(get_size(n.name))
        push!(all_sizes,get_size(n.name))
    end
    tot_size = sum(all_sizes)

    label_color_dict = Dict()
    for n in getleaflist(newt)
        if get_size(n.name) > tot_size/10
            label_color_dict[n] = "red"
        else
            label_color_dict[n] = "black"
        end
    end

    img = tree_draw(newt,dot_color_default = "black",draw_labels = true,label_color_dict = label_color_dict,
        line_width = 0.3, dot_size_dict = dot_size_dict,
        max_dot_size = 0.03, dot_color_dict = label_color_dict,
        dot_opacity = 0.85,
    canvas_height = length(getleaflist(newt))*0.2cm, canvas_width = 10cm
    )
    img |> SVG(filepath*".svg", 10cm, length(getleaflist(newt))*0.2cm)
end

const AA_colors = [
    'Q' => RGB(1.0,0.0,0.8),
    'E' => RGB(1.0,0.0,0.4),
    'D' => RGB(1.0,0.0,0.0),
    'S' => RGB(1.0,0.2,0.0),
    'T' => RGB(1.0,0.4,0.0),
    'G' => RGB(1.0,0.6,0.0),
    'P' => RGB(1.0,0.8,0.0),
    'C' => RGB(1.0,1.0,0.0),
    'A' => RGB(0.8,1.0,0.0),
    'V' => RGB(0.6,1.0,0.0),
    'I' => RGB(0.4,1.0,0.0),
    'L' => RGB(0.2,1.0,0.0),
    'M' => RGB(0.0,1.0,0.0),
    'F' => RGB(0.0,1.0,0.4),
    'Y' => RGB(0.0,1.0,0.8),
    'W' => RGB(0.0,0.8,1.0),
    'H' => RGB(0.0,0.4,1.0),
    'R' => RGB(0.0,0.0,1.0),
    'K' => RGB(0.4,0.0,1.0),
    'N' => RGB(0.8,0.0,1.0),
    '-' => "grey60"
];

const NT_colors = [
      'A' => "salmon",
      'G' => "gold",
      'T' => "lightgreen",
      'C' => "lightblue",
      '-' => "grey60",
]

function highlighter_figure(fasta_collection; out_path = "figure.png")
    #read in and collapse
    seqnames, ali_seqs = read_fasta_with_names(fasta_collection);
    collapsed_seqs, collapsed_sizes, collapsed_names = variant_collapse(ali_seqs; prefix = "v")

    if length(collapsed_seqs) == 1
        @warn "All sequences in $(split(basename(fasta_collection),".fast")[1]) are identical, so no collapsed tree can be created"
        SVG(out_path, 20cm, 20cm) #blank SVG
    else
        treestring = fasttree_nuc(collapsed_seqs, collapsed_names; quiet = true)[1];
        newt = MolecularEvolution.gettreefromnewick(treestring, MolecularEvolution.GeneralFelNode);
        root = [n for n in MolecularEvolution.getleaflist(newt) if split(n.name,'_')[1] == "v1"][1]
        rerooted = MolecularEvolution.ladderize(MolecularEvolution.reroot(root))

        dot_size_dict = Dict()
        all_sizes = []
        for n in MolecularEvolution.getleaflist(rerooted)
            dot_size_dict[n] = sqrt(get_size(n.name))
            push!(all_sizes,get_size(n.name))
        end
        tot_size = sum(all_sizes)

        color_dict = Dict()
        for n in MolecularEvolution.getleaflist(rerooted)
            if get_size(n.name) > tot_size/10
                color_dict[n] = "red"
            else
                color_dict[n] = "black"
            end
        end

        img = MolecularEvolution.highlighter_tree_draw(rerooted, uppercase.(collapsed_seqs), collapsed_names, uppercase(collapsed_seqs[1]);
            legend_colors = NT_colors, legend_padding = 1cm,
            tree_args = [:line_width => 0.5mm, :font_size => 1, :dot_size_dict => dot_size_dict, :label_color_dict => color_dict,
                :dot_color_dict => color_dict, :max_dot_size => 0.1, :dot_opacity => 0.6, :name_opacity => 0.8])
        compose(context(0.05,0.01,0.95,0.99), img) |> SVG(out_path, 20cm, 20cm)
    end
end

#-
function gettreefromnewick(str, T::DataType; tagged = false, disable_binarize = false)
    currnode = T()
    i = 1

    str = collect(str)

    function try_apply_char_arr(node, char_arr)
        if node.name == ""
            #println(char_arr)
            node.name = join(char_arr)
            #println("naming a node: ", join(char_arr))
        else
            #println("yo!")
        end
        empty!(char_arr)
    end
    char_arr = []
    tag_dict = Dict{T,Vector{String}}()
    while i <= length(str)
        c = str[i]
        
        if c == '{' && tagged
            init_loc = (i += 1)
            while i <= length(str)
                #println(i)
                if (str[i] == '}')
                    break
                end
                i += 1
            end
            tag_dict[currnode] = string.(split(join(str[init_loc:i-1]), ","))
            i += 1
        elseif c == '('
            #println("making new node")
            try_apply_char_arr(currnode, char_arr)
            newnode = T()
            addchild(currnode, newnode)
            currnode = newnode
            i += 1
        elseif c == ')'
            try_apply_char_arr(currnode, char_arr)
            currnode = currnode.parent
            i += 1

        elseif c == ':'
            try_apply_char_arr(currnode, char_arr)
            valid_nums = ['1','2','3','4','5','6','7','8','9','0','-','.', 'E', 'e']
            init_loc = (i += 1)
            while i <= length(str)
                #println(i)
                if !(str[i] in valid_nums)
                    break
                end
                i += 1
            end
            currnode.branchlength = parse(Float64, join(str[init_loc:i-1]))

        elseif c == ','
            #println("making new node")
            try_apply_char_arr(currnode, char_arr)
            newnode = T()
            addchild(currnode.parent, newnode)
            currnode = newnode
            i += 1
        elseif c == ';'
            try_apply_char_arr(currnode, char_arr)
            return (tagged ? (currnode, tag_dict) : currnode)
        else
            push!(char_arr, c)
            #println(char_arr)
            i += 1
        end
    end

    binarize!(currnode)
    
    return (tagged ? (currnode, tag_dict) : currnode)
end
-#
