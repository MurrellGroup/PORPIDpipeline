
"""
    function H704_init_template_proc(fasta_collection, panel_file,
      mds_plot_file, apobec_file, realigned_tre_plot_file, realigned_file;
      agreement_thresh=0.7,panel_thresh=50, af_thresh=0.15)
Given a fasta collection, filter, produce and write mafft alignments,
MDS plots with APOBEC model, and save MolEv phylogeny.
"""
function H704_init_template_proc(fasta_collection, panel_file,
      mds_plot_file, apobec_file, realigned_tre_plot_file, realigned_file;
      agreement_thresh=0.7, panel_thresh=50, af_thresh=0.15) #, af_cutoff=1)
    reject_df = DataFrame(reject_reason=String[],threshold=Float64[],count=Int64[])
    f = fasta_collection
    #align consensus
    seqs,seqnames,sizes,agreement_scores = read_HVTN(f)
    annot_names = seqnames .* " num_CCS=" .* string.(sizes) .* " min_agreement=" .* string.(agreement_scores)
    reject_seqs=[]
    reject_names=[]
    
    # filter reads that do not make minimum agreement
    push!(reject_df,["minimum_agreement < threshold",agreement_thresh,
        sum(agreement_scores .< agreement_thresh)])
    
    reject_seqs = vcat(reject_seqs, seqs[agreement_scores .< agreement_thresh])
    reject_names = vcat(reject_names, annot_names[agreement_scores .< agreement_thresh] .* " no_min_agreement" )
    
    #filter, minimum_agreement >= agreement_thresh (default 0.7)
    seqs = seqs[agreement_scores .>= agreement_thresh]
    annot_names = annot_names[agreement_scores .>=agreement_thresh]
    sizes = sizes[agreement_scores .>= agreement_thresh]
    agreement_scores = agreement_scores[agreement_scores .>= agreement_thresh]
    
    # calculate af_cutoff after discarding seq that dont meet minimum agreement
    af_cutoff = artefact_cutoff(sizes,af_thresh)
    
    # filter artefacts
    push!(reject_df,["ccs_count < artefact cutoff",af_cutoff,sum(sizes .< af_cutoff)])
    reject_seqs=vcat(reject_seqs, seqs[sizes .< af_cutoff])
    reject_names=vcat(reject_names, annot_names[sizes .< af_cutoff].*" possible_artefact_at_$(af_cutoff)_cutoff" )
    
    seqs = seqs[sizes .>= af_cutoff]
    annot_names = annot_names[sizes .>= af_cutoff]
    agreement_scores = agreement_scores[sizes .>= af_cutoff]
    sizes = sizes[sizes .>= af_cutoff]
    
    ali_seqs = mafft_align(seqs);
    
    #clean seqs
    ref_names, ref_seqs = read_fasta(panel_file)
    ref_profile = seqs2profile(ref_seqs)
    extracted_seqs, scores = extract_and_score_misalignments(ali_seqs, ref_profile)
    
    #realign
    realigned_seqs = mafft_align(extracted_seqs[scores .< panel_thresh]);
    realigned_names = annot_names[scores .< panel_thresh]
    realigned_sizes = sizes[scores .< panel_thresh]
    cons = ""
    if length(realigned_seqs) == 0
        @error "No sequences from $(split(basename(f),".fast")[1]) pass the panel threshold. Check sequences in contam_passed folder and alter panel_thresh in Snakefile if needed"
    else
        cons = gap_preserving_consensus(realigned_seqs)
    end
    
    #probabilistic APOBEC model
    model_results = APOBEC.(cons,realigned_seqs)
    post_APOBEC = [(el[2]) for el in model_results];
    apobec_df = modelresults2table(model_results, realigned_names)
    CSV.write(apobec_file,apobec_df)

    #MDS with p(APOBEC|mutations)
    tight_layout()
    fig = ali2MDS(realigned_seqs,
            realigned_sizes,
            post_APOBEC)
    title("P(APOBEC|mutations): $(split(basename(f),'.')[1])")
    xlab, ylab = xlabel("MDS 1"), ylabel("MDS 2")

    #write alignment (IO step)
    write_fasta(realigned_file,
        realigned_seqs,
        names = realigned_names)

    push!(reject_df,["dist_from_panel >= threshold",panel_thresh,sum(scores .>= panel_thresh)])
    reject_seqs=vcat(reject_seqs,extracted_seqs[scores .>= panel_thresh])
    reject_names=vcat(reject_names,annot_names[scores .>= panel_thresh] .* " panel_mismatch")
    write_fasta(realigned_file*".rejected.fasta",reject_seqs,names=reject_names)
        # extracted_seqs[scores .>= panel_thresh]; # removed degap.()
        # names = annot_names[scores .>= panel_thresh])
    # panel_reject_df = DataFrame(Theshold_Percentage=[panel_thresh],
    #                            Greater_Than_Theshold=[sum(scores .>= panel_thresh)])
    CSV.write(realigned_file*".rejected.csv",reject_df)

    #save MDS plot (IO step)
    fig.savefig(mds_plot_file;
        transparent = true,
        dpi = 200,
        bbox_inches = "tight")

    #highlighter plot and tree
    highlighter_figure(realigned_file; out_path = realigned_tre_plot_file)
    return realigned_seqs, realigned_names, af_cutoff
end

"""
  function H704_ENV_proc(fasta_path,ali_seqs,seqnames;
    env_profile=env_profile, re_profile=re_profile)
      ---- Not in use ------
"""
function H704_ENV_proc(fasta_path,ali_seqs,seqnames;
      env_profile=env_profile, re_profile=re_profile)
    f = fasta_path
    sample_id = join(split(basename(f),'_')[1:3],'_');

    env_seqs,(env_start,env_end) = extract_region(ali_seqs,env_profile);
    write_fasta(snakemake.output[4]*".env.fa",env_seqs,names = seqnames)

    re_seqs,(re_start,re_end) = extract_region(ali_seqs, re_profile);
    write_fasta(snakemake.output[4]*".re.fa", re_seqs; names = seqnames);

    #First function checker for alignment based start, but first occuring stop codon, no matter where it happens.
    from_start_of_env_AA = robust_translate.(a[env_start:end] for a in ali_seqs)
    function_checked = env_function_check_2.(from_start_of_env_AA)
    reason_names = seqnames .* " " .* [i[3] for i in function_checked]
    translated = [i[2] for i in function_checked]
    keeps = [i[1] for i in function_checked]
    AA_keeps_ali = mafft_stop_preserving_align(translated[keeps])
    write_fasta(snakemake.output[4]*".envORF.AA.passed.aligned.fa",AA_keeps_ali,names = seqnames[keeps])
    orf_synth_names, orf_synth_seqs = pick_RE_synth_variants_by_thresh(ali_seqs[keeps],env_start,length(ali_seqs[1]),re_start,length(ali_seqs[1]);thresh=0.095)

    #find native stop
    stop_pos = find_first_stop_pos_with_gaps.(s[env_start:end] for s in orf_synth_seqs);
    trimmed_orf_synth_seqs = []
    for (i,seq) in enumerate(orf_synth_seqs)
        trim = seq[re_start:env_start + stop_pos[i] - 1]
        push!(trimmed_orf_synth_seqs, trim)
    end
    write_fasta(snakemake.output[4]*".envORF.RE.synth.fa",degap.(trimmed_orf_synth_seqs);names="$(sample_id)_".*orf_synth_names)
    write_fasta(snakemake.output[4]*".envORF.AA.failed.fa",from_start_of_env_AA[.!keeps],names = reason_names[.!keeps])
    write_fasta(snakemake.output[4]*".envORF.nuc.failed.fa",[a[env_start:end] for a in ali_seqs][.!keeps],names = reason_names[.!keeps])

    #Second function checker based on reference-based start and end coords.
    AA_seqs = robust_translate.(env_seqs);
    function_checked = env_function_check.(AA_seqs)
    keeps = [i[1] for i in function_checked]
    translated = [i[2] for i in function_checked]
    AA_keeps_ali = mafft_stop_preserving_align(translated[keeps])
    synth_names, synth_seqs = pick_RE_synth_variants_by_thresh(ali_seqs[keeps],env_start,env_end,re_start,re_end;thresh=0.095)
    write_fasta(snakemake.output[4]*".env.RE.synth.fa",degap.(synth_seqs);names="$(sample_id)_".*synth_names)

    write_fasta(snakemake.output[4]*".env.AA.passed.aligned.fa",AA_keeps_ali,names = seqnames[keeps])
    ali_env_keeps = mafft_align(env_seqs[keeps])
    write_fasta(snakemake.output[4]*".env.nuc.passed.aligned.fa",ali_env_keeps,names = seqnames[keeps])
    #Filtered Trimmed Nuc Variant collapse and tree drawing
    Vseqs,Vsizes,Vnames = variant_collapse(ali_env_keeps,prefix = "FTNV") #Filtered Trimmed Nuc Variant
    size_dict = Dict(zip(Vnames,Vsizes))
    treestring = fasttree_nuc(Vseqs,Vnames; quiet = true)[1]
    draw_variant_tree(treestring, leaf_rename, snakemake.output[4]*".filtered.env.tre")

    reason_names = seqnames .* " " .* [i[3] for i in function_checked]
    write_fasta(snakemake.output[4]*".env.nuc.failed.aligned.fa",env_seqs[.!keeps],names = reason_names[.!keeps])
    write_fasta(snakemake.output[4]*".env.AA.failed.fa",AA_seqs[.!keeps],names = reason_names[.!keeps])
end

