using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, StatsBase
using BioSequences, FASTX
using DataFrames, CSV
using CodecZlib: GzipDecompressorStream
using CodecZlib: GzipCompressorStream

println("using Julia version: $(VERSION)")

t1 = time()

SAMPLE_CONFIGS = snakemake.params["config"]
mkdir(snakemake.output[1])

verbose = false
if snakemake.params["verbose"] == "true"
    verbose = true
end

f_kwargs = [
    :demux_dir => snakemake.output[1],
    :samples => SAMPLE_CONFIGS,
    :verbose => verbose,
    :error_rate => snakemake.params["error_rate"],
    :min_length => snakemake.params["min_length"],
    :max_length => snakemake.params["max_length"],
    :label_prefix => "seq",
    :error_out => true
    ]

println("performing chunked quality filtering and demux on $(snakemake.input[1])")
chunk_size = snakemake.params["chunk_size"]

reads = chunked_filter_apply(snakemake.input[1], snakemake.output[1], chunked_quality_demux;
    chunk_size=chunk_size, f_kwargs)

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
      close(stream)
    end
    sample_name = replace( replace( basename(path), ".gz" => "" ), ".fastq" => "" )
    println(sample_name," => ",count)
    if occursin("REJECT",path)
        push!(df_reject,[sample_name,count])
    else
        push!(df_demux,[sample_name,count])
    end
end
println("-------------------------------------")
println("total demuxed => $(no_assigned)")
println("total rejected => $(no_rejected)")

# now do the down sampling
df_demux_sampled = DataFrame(Sample = [], Count = [], Sampled = [])
max_reads = snakemake.params["max_reads"]
if max_reads < 1
    max_reads = 10000000
end
println("downsampling to about $(max_reads) reads")
global no_lost=0
global no_retained=0
for path in filepaths
    out_path = path[1:end-9]*"_sampled.fastq.gz"
    if ! occursin("REJECT",path)
        if endswith(path, ".gz")
            stream = FASTQ.Reader(GzipDecompressorStream(open(path)))
            out_stream = FASTQ.Writer(GzipCompressorStream(open(out_path,"w")))
        else
            stream = FASTQ.Reader(open(path))
            out_stream = FASTQ.Writer(open(out_path,"w"))
        end
        sample_name = replace( replace( basename(path), ".gz" => "" ), ".fastq" => "" )
        count = 0
        sampled=0
        ac = df_demux[df_demux.Sample .== sample_name, :Count][1]
        sp = 1.0
        ac > 0 ? sp = min(1.0, max_reads / ac) : sp = 1.0
        if sp < 1.0
            for record in stream
                count += 1
                if rand() < sp
                    write(out_stream, record)
                    sampled += 1
                    global no_retained += 1
                else
                    global no_lost += 1
                end
            end
            close(stream)
            close(out_stream)
            println("$(sample_name), sampled $(sampled) reads")
            rm(path)
            mv(out_path,path)
            push!(df_demux_sampled,[sample_name,count,sampled])
        else
            close(stream)
            close(out_stream)
            println("$(sample_name), retained all $(ac) reads")
            global no_retained += ac
            rm(out_path)
            push!(df_demux_sampled,[sample_name,ac,ac])
        end
    end
end

CSV.write("$(snakemake.output[3])", df_demux_sampled)
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
push!(df_qual,["Number lost to downsampling",no_lost])
push!(df_qual,["Number retained after downsampling",no_retained])
CSV.write("$(snakemake.output[2])", df_qual)


t2 = time()
println("Quality filtering and demultiplexing took $(t2 - t1) seconds.")
