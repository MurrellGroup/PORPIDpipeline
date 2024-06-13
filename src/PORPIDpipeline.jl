module PORPIDpipeline

export # apobec_model
    APOBEC
export # functions
    mafft,
    mafft_align,
    family_size_umi_len_stripplot,
    family_size_stripplot,
    highlighter_figure,
    gettreefromnewick,
    di_nuc_freqs,
    artefact_cutoff
# export # molev_functions
#     highlighter_figure
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
    chunked_quality_demux
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
export # blast_functions
    Hamming,
    PairWise,
    get_blast_results

using BioSequences, FASTX

include("apobec_model.jl")
include("functions.jl")
# include("molev_functions.jl")
include("porpid_analysis_methods.jl")
include("postproc_functions.jl")
include("demux_functions.jl")
include("contam-filter_functions.jl")
include("porpid_functions.jl")
include("blast_functions.jl")

end # module
