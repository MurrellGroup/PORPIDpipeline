import subprocess, sys
configfile: "config.yaml"
DATASETS = [d for d in config for s in config[d]]
SAMPLES = [s for d in config for s in config[d]]
VERSION = "1.10.5"
COMMIT = subprocess.check_output(['git', 'rev-parse', '--verify', 'HEAD']).strip().decode()

sys.stderr.write("Running PorpidPostproc\n")
sys.stderr.write("Version: {0}\n".format(VERSION))
sys.stderr.write("Commit ID: {0}\n".format(COMMIT))

def contam_input(wildcards):
    SAMPLES = [s for s in config[wildcards.dataset]]
    return expand("porpid/{dataset}/consensus/{sample}.fasta",
        dataset = wildcards.dataset,
        sample = SAMPLES
    )

def bottle1_input(wildcards):
    SAMPLES = [s for s in config[wildcards.dataset]]
    return expand("porpid/{dataset}/tags/{sample}.csv",
        dataset = wildcards.dataset,
        sample = SAMPLES
    )
    
def bottle2_input(wildcards):
    SAMPLES = [s for s in config[wildcards.dataset]]
    return expand("postproc/{dataset}/{sample}/{sample}.fasta",
        dataset = wildcards.dataset,
        sample = SAMPLES
    )

# PorpidPostproc parameters
# demux
chunk_size = 100000   # default 10000
error_rate = 0.01    # default 0.01
min_length = 2100    # default 2100
max_length = 4300    # default 4300
#porpid
fs_thresh = 5        # default 1 (must be 1 for artefact filter to work)
lda_thresh = 0.995   # default 0.995
#contam
cluster_thresh = 0.015   # default 0.015
proportion_thresh = 0.2  # default 0.2
dist_thresh = 0.015      # default 0.015
# change to "off" to switch off contam filter
contam_toggle = "on"   # default "on" 
#postproc
af_thresh = 0.0  # default 0.15 (drops smallest 15% of CCS reads)
agreement_thresh = 0.7   # default 0.7
panel_thresh = 50        # default 50
#tar
# change this to "false" if you dont want 
# degapped postproc sequences
degap = "true"           # default "true"

rule all:
    input:
        expand("porpid/{dataset}-porpid.tar.gz",
            dataset = DATASETS),
        expand("postproc/{dataset}-postproc.tar.gz",
            dataset = DATASETS)

rule demux:
    input:
        "raw-reads/{dataset}.fastq.gz"
    output:
        directory("porpid/{dataset}/demux"),
        "porpid/{dataset}/quality_report.csv",
        "porpid/{dataset}/demux_report.csv",
        "porpid/{dataset}/reject_report.csv"
    params:
        chunk_size = chunk_size,
        error_rate = error_rate,
        min_length = min_length,
        max_length = max_length,
        config = lambda wc: config[wc.dataset]
    script:
        "scripts/demux.jl"

rule porpid:
    input:
        "porpid/{dataset}/demux"
    output:
        directory("porpid/{dataset}/porpid/{sample}.fastq"),
        "porpid/{dataset}/tags/{sample}.csv"
        # directory("porpid/{dataset}/porpid/{sample}.fastq/{sample}")
    params:
        config = lambda wc: config[wc.dataset][wc.sample],
        fs_thresh = fs_thresh,
        lda_thresh= lda_thresh
    script:
        "scripts/porpid.jl"
        
rule porpid_wait:
    input:
        files = bottle1_input
    output:
        "porpid/{dataset}/bottle1_report.csv"
    script:
        "scripts/bottle.jl"

rule consensus:
    input:
        "porpid/{dataset}/porpid/{sample}.fastq",
        "porpid/{dataset}/tags/{sample}.csv",
        "porpid/{dataset}/bottle1_report.csv"
    output:
        "porpid/{dataset}/consensus/{sample}.fasta",
        "porpid/{dataset}/tags_filtered/{sample}.csv"
    params:
        config = lambda wc: config[wc.dataset][wc.sample],
        af_thresh = af_thresh,
        agreement_thresh = agreement_thresh
    script:
        "scripts/consensus.jl"

rule contam:
    input:
        files = contam_input,
        panel = "panels/contam_panel.fasta"
    output:
        directory("porpid/{dataset}/contam_passed"),
        directory("porpid/{dataset}/contam_failed"),
        "porpid/{dataset}/contam_report.csv",
        "porpid/{dataset}/contam_suspect.csv"
    params:
        proportion_thresh =proportion_thresh,
        cluster_thresh = cluster_thresh,
        dist_thresh = dist_thresh,
        contam_toggle = contam_toggle
    script:
        "scripts/contam.jl"

rule postproc:
    input:
        "porpid/{dataset}/contam_passed",
        "porpid/{dataset}/tags_filtered/{sample}.csv",
        "porpid/{dataset}/porpid/{sample}.fastq"
    output:
        report("postproc/{dataset}/{sample}/{sample}.fasta.mds.png", category = "postproc", caption = "report-rst/mds.rst"),
        "postproc/{dataset}/{sample}/{sample}.fasta.apobec.csv",
        report("postproc/{dataset}/{sample}/{sample}.fasta.tre.svg", category = "postproc", caption = "report-rst/highlighter.rst"),
        "postproc/{dataset}/{sample}/{sample}.fasta",
        report("postproc/{dataset}/{sample}/{sample}_qc_bins.png", category = "postproc", caption = "report-rst/bins.rst"),
        report("postproc/{dataset}/{sample}/{sample}_artefacts.png", category = "postproc", caption = "report-rst/artefacts.rst"),
        "postproc/{dataset}/{sample}/{sample}_qc_bins.csv",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.fasta",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.csv",
        report("postproc/{dataset}/{sample}/{sample}_di_nuc_freq.png", category = "postproc", caption = "report-rst/di_nuc_freq.rst")
    params:
        panel = lambda wc: config[wc.dataset][wc.sample]["panel"],
        af_thresh = af_thresh,
        fs_thresh = fs_thresh,
        agreement_thresh = agreement_thresh,
        panel_thresh = panel_thresh
    script:
        "scripts/postproc.jl"
        
rule postproc_wait:
    input:
        files = bottle2_input
    output:
        "postproc/{dataset}/bottle2_report.csv"
    script:
        "scripts/bottle.jl"

rule report:
    input:
        "postproc/{dataset}/{sample}/{sample}_qc_bins.png",
        "postproc/{dataset}/{sample}/{sample}_artefacts.png",
        "postproc/{dataset}/{sample}/{sample}_qc_bins.csv",
        "postproc/{dataset}/{sample}/{sample}.fasta.mds.png",
        "postproc/{dataset}/{sample}/{sample}.fasta.tre.svg",
        "postproc/{dataset}/{sample}/{sample}.fasta",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.fasta",
        "postproc/{dataset}/{sample}/{sample}.fasta.rejected.csv",
        "postproc/{dataset}/{sample}/{sample}_di_nuc_freq.png",
        "postproc/{dataset}/bottle2_report.csv"
    params:
        VERSION = VERSION,
        COMMIT = COMMIT
    output:
        "postproc/{dataset}/{sample}/{sample}-report.html"
    script:
        "scripts/report.jl"
        
rule index:
    input:
        expand("postproc/{dataset}/{sample}/{sample}-report.html", zip, dataset = DATASETS, sample = SAMPLES),
        expand("postproc/{dataset}/{sample}/{sample}_qc_bins.csv", zip, dataset = DATASETS, sample = SAMPLES)
    params:
        VERSION = VERSION,
        COMMIT = COMMIT,
        SAMPLES = SAMPLES,
        chunk_size = chunk_size,
        error_rate = error_rate,
        min_length = min_length,
        max_length = max_length,
        proportion_thresh =proportion_thresh,
        cluster_thresh = cluster_thresh,
        dist_thresh = dist_thresh,
        fs_thresh = fs_thresh, 
        af_thresh = af_thresh,
        lda_thresh = lda_thresh,
        agreement_thresh = agreement_thresh,
        panel_thresh = panel_thresh,
        contam_toggle = contam_toggle
    output:
        "postproc/{dataset}/{dataset}-index.html",
        "postproc/{dataset}/{dataset}-sequence_report.csv"
    script:
        "scripts/index.jl"
        
rule tar:
    input:
        "postproc/{dataset}/{dataset}-index.html"
    output:
        "porpid/{dataset}-porpid.tar.gz",
        "postproc/{dataset}-postproc.tar.gz"
    params:
        degap = degap,
        datasets = DATASETS,
        samples = SAMPLES
    script:
        "scripts/tar.jl"
