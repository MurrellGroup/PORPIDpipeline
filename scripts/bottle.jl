using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using DataFrames,CSV, DataFramesMeta

files = snakemake.input["files"]

df = DataFrame(completed=files)
CSV.write(snakemake.output[1],df)
println("bottleneck passed after $(length(files)) processes completed...")
