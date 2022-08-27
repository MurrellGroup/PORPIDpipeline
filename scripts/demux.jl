ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, NextGenSeqUtils, StatsBase
using BioSequences, FASTX
using DataFrames, CSV
using CodecZlib: GzipDecompressorStream
# include("../../src/postproc_functions.jl")

println("using Julia version: $(VERSION)")

t1 = time()

SAMPLE_CONFIGS = snakemake.params["config"]
mkdir(snakemake.output[1])

f_kwargs = [
    :demux_dir => snakemake.output[1],
    :samples => SAMPLE_CONFIGS,
    :verbose => false,
    :error_rate => snakemake.params["error_rate"],
    :min_length => snakemake.params["min_length"],
    :max_length => snakemake.params["max_length"],
    :label_prefix => "seq",
    :error_out => true
    ]

println("performing chunked quality filtering and demux on $(snakemake.input[1])")
chunk_size = snakemake.params["chunk_size"]

reads = chunked_filter_apply(snakemake.input[1], snakemake.output[1], chunked_quality_demux;
    chunk_size=chunk_size, f_kwargs, verbose = false)

# create empty fastq files for those samples with no reads
# this does not work, too much trouble with downstream scripts
# for sample in snakemake.params["SAMPLES"]
#     println(sample)
#     run(`touch "$(snakemake.output[1])/$(sample).fastq"`)
# end

# now do some accounting
total_reads = reads[1]
quality_reads = reads[2]
bad_reads = reads[3]
short_reads = reads[4]
long_reads = reads[5]
demuxed_reads = reads[6]

# report on total sequences for each sample
println()
println("total reads => $(total_reads)")
println("quality reads => $(quality_reads)")
println("bad reads => $(bad_reads)")
println("short reads => $(short_reads)")
println("long reads => $(long_reads)")
println("demuxed reads => $(demuxed_reads)")
println("-------------------------------------")
no_assigned = 0
no_rejected = 0
filepaths = readdir(snakemake.output[1],join=true)
df_demux = DataFrame(Sample = [], Count = [])
df_reject = DataFrame(RejectFile = [], Count = [])
for path in filepaths
    if endswith(path, ".gz")
        stream = FASTQ.Reader(GzipDecompressorStream(open(path)))
    else
        stream = FASTQ.Reader(open(path))
    end
    # stream = open(FASTQ.Reader, path)
    count = 0
    if filesize(path) > 0
      for record in stream
        count += 1
        if occursin("REJECT",path)
            global no_rejected += 1
        else
            global no_assigned += 1
        end
      end
    end
    # sample_name = replace(basename(path), ".fastq" => "")
    sample_name = basename(path)
    println(sample_name," => ",count)
    if occursin("REJECT",path)
        push!(df_reject,[sample_name,count])
    else
        push!(df_demux,[replace(sample_name,".fastq" => ""),count])
    end
end
println("-------------------------------------")
println("total demuxed => $(no_assigned)")
println("total rejected => $(no_rejected)")
CSV.write("$(snakemake.output[3])", df_demux)
CSV.write("$(snakemake.output[4])", df_reject)

# save a quality_filter report
# no_of_fails = no_of_reads - no_assigned # fix this, get chunked_filter_apply to return no_of_fails
df_qual = DataFrame(Description = [], Count = [])
push!(df_qual,["Initial number of raw reads",total_reads])
push!(df_qual,["Number of quality reads",quality_reads])
push!(df_qual,["Number of bad reads",bad_reads])
push!(df_qual,["Number of short reads",short_reads])
push!(df_qual,["Number of long reads",long_reads])
push!(df_qual,["Number assigned by demux",no_assigned])
CSV.write("$(snakemake.output[2])", df_qual)


t2 = time()
println("Quality filtering and demultiplexing took $(t2 - t1) seconds.")
