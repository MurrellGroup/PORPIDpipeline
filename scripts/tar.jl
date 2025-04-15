using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

using PORPIDpipeline, FASTX

# zip porpid and postproc directories for easy download

dataset = snakemake.wildcards["dataset"]
porpid_dir = "porpid/$(dataset)"
postproc_dir = "postproc/$(dataset)"

degap_flag = snakemake.params["degap"]
println( "degap flag = $(degap_flag)")
samples = snakemake.params["samples"]
if degap_flag == "true"
    t1s = time()
    println("generating degapped postproc sequences...")
    run(`mkdir -p $(postproc_dir)/degapped_fasta/`)
    for sample in samples
        infile = "$(postproc_dir)/$(sample)/$(sample).fasta"
        outfile = "$(postproc_dir)/degapped_fasta/$(sample)-degapped.fasta"
        print(".")
        names, descripts, seqs = read_fasta_with_names_and_descriptions(infile)
        write_fasta(outfile, degap.(seqs), names = descripts)
    end
    println()
    t1f = time()
    println("degapping took $(t1f - t1s) seconds.")
end

collapse_flag = snakemake.params["collapse"]
println( "collapse flag = $(collapse_flag)")
samples = snakemake.params["samples"]
if collapse_flag == "true"
    t2s = time()
    println("generating collapsed postproc sequences...")
    run(`mkdir -p $(postproc_dir)/collapsed_fasta/`)
    for sample in samples
        infile = "$(postproc_dir)/$(sample)/$(sample).fasta"
        outfile = "$(postproc_dir)/collapsed_fasta/$(sample)-collapsed.fasta"
        print(".")
        names, descripts, seqs = read_fasta_with_names_and_descriptions(infile)
        col_seqs, col_sizes, col_names = variant_collapse(seqs,
            prefix = "$(sample)_v")
        write_fasta(outfile, col_seqs, names = col_names)
    end
    println()
    t2f = time()
    println("collapsing took $(t2f - t2s) seconds.")
end

# and now tar and zip both porpid and postproc directories
# first copy some reports from porpid to postproc
run(`cp $(porpid_dir)/contam_report.csv $(postproc_dir)/contam_report.csv`)

# now archive the postproc_dir
# first rename postproc_dir, tar and zip and then rename back
t3s = time()
println("archiving postproc ...")
run(`mv $(postproc_dir) $(postproc_dir)-postproc`)
run(`tar -C postproc -czf postproc/$(dataset)-postproc.tar.gz $(dataset)-postproc`)
run(`mv $(postproc_dir)-postproc $(postproc_dir)`)
t3f = time()
println("postproc archiving took $(t3f - t3s) seconds.")

porpid_archive_flag = snakemake.params["porpid_archive"]
if porpid_archive_flag == "full"
    # now rename porpid_dir, tar and zip and rename back
    println("archiving FULL porpid directory ...")
    t4s = time()
    run(`mv $(porpid_dir) $(porpid_dir)-porpid`)
    run(`tar -C porpid -czf porpid/$(dataset)-porpid.tar.gz $(dataset)-porpid`)
    run(`mv $(porpid_dir)-porpid $(porpid_dir)`)
    t4f = time()
    println("full porpid archiving took $(t4f - t4s) seconds.")
else
    # now rename porpid_dir, tar and zip and rename back
    println("archiving PARTIAL porpid directory ...")
    t5s = time()
    run(`mkdir -p $(porpid_dir)-porpid`)
    run(`cp $(porpid_dir)/demux_report.csv $(porpid_dir)-porpid`)
    run(`cp $(porpid_dir)/quality_report.csv $(porpid_dir)-porpid`)
    run(`cp $(porpid_dir)/reject_report.csv $(porpid_dir)-porpid`)
    run(`cp $(porpid_dir)/contam_report.csv $(porpid_dir)-porpid`)
    run(`cp $(porpid_dir)/contam_suspect.csv $(porpid_dir)-porpid`)
    run(`cp -r $(porpid_dir)/demux $(porpid_dir)-porpid`)
    run(`cp -r $(porpid_dir)/tags $(porpid_dir)-porpid`)
    run(`cp -r $(porpid_dir)/tags_filtered $(porpid_dir)-porpid`)
    run(`cp -r $(porpid_dir)/consensus $(porpid_dir)-porpid`)
    run(`cp -r $(porpid_dir)/contam_passed $(porpid_dir)-porpid`)
    run(`cp -r $(porpid_dir)/contam_failed $(porpid_dir)-porpid`)
    run(`tar -C porpid -czf porpid/$(dataset)-porpid.tar.gz $(dataset)-porpid`)
    run(`rm -rf $(porpid_dir)-porpid`)
    t5f = time()
    println("partial porpid archiving took $(t5f - t5s) seconds.")
end

