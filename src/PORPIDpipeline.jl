module PORPIDpipeline

export # apobec_model
    APOBEC
export # functions
    mafft,
    mafft_align,
    family_size_umi_len_stripplot,
    family_size_stripplot,
    variant_collapse,
    highlighter_figure,
    gettreefromnewick,
    di_nuc_freqs,
    artefact_cutoff
export # porpid-analysis-methods
    generateConsensusFromDir
export # postproc_functions
    H704_init_template_proc
export # demux_functions
    unique_not_substr,
    longest_conserved_5p,
    iterative_primer_match,
    sliding_demux_dic,
    chunked_filter_apply,
    chunked_quality_demux,
    fast_primer_match,
    fast_primer_pair_match,
    demux_dict,
    primer_trim,
    primer_trim_reverse,
    double_primer_trim,
    primer_peek,
    chunked_fastq_filter_demux

export # contam-filter_functions
    IUPACbool,
    resolve_base,
    resolve_seq,
    db_cluster_seqs,
    db_seqs,
    contam_check,
    read_fasta_with_descriptors_in_names
export # porpid_functions
    filterCCSFamilies,
    porpid_write_to_file,
    porpid_write_to_dictionary,
    porpid_write_to_file_count_to_dict
export # pipeline_utils (lifted from NextGenSeqUtils)
    read_fastq,
    write_fastq,
    read_fasta_records,
    read_fasta,
    write_fasta,
    read_fasta_with_names_and_descriptions,
    degap,
    reverse_complement,
    freq,
    maxfreq,
    kmer_count,
    corrected_kmer_dist,
    PATHS,
          # export all of align.jl from NextGenSeqUtils
    usearch_filter,
    usearch_trim_fastq_with_phreds,

    nw_align,
    banded_nw_align,
    triplet_nw_align,
    local_align,
    kmer_seeded_align,
    triplet_kmer_seeded_align,
    loc_kmer_seeded_align,
    local_kmer_seeded_align,
    kmer_seeded_edit_dist,
    resolve_alignments,
    align_reference_frames,
    local_edit_dist,
    merge_alignments,
    seqs2profile,
    profile_affine_align,
    gap_elem,
    profile_cost,
    affine_nw_align


using BioSequences, FASTX

include("apobec_model.jl")
include("functions.jl")
include("porpid_analysis_methods.jl")
include("postproc_functions.jl")
include("demux_functions.jl")
include("contam-filter_functions.jl")
include("porpid_functions.jl")
include("pipeline_utils.jl")

end # module
