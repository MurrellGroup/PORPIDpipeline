ENV["MPLBACKEND"] = "Agg"
using Base64, CSV, DataFrames, NextGenSeqUtils, FASTX, WebBlast

function Hamming(s1, s2)
    if length(s1) != length(s2)
        @warn "cant compute Hamming dist on strings of unequal length"
        return max(length(s1),length(s2))
    end
    s1v=split(s1,"")
    s2v=split(s2,"")
    sgv=["-" for i in 1:length(s1)]
    return sum( (s1v .!= s2v) .& (s1v .!= sgv) .& (s2v .!= sgv) )
end

function PairWise(f,xs,ys)
    return [f(x,y) for x in xs, y in ys]
end

function get_blast_results(files; thresh_hold=0.05, max_clades=5, max_waits=30, verbose=true)
    queries = DataFrame(sample=String[], count=Int[], query=String[])
    for file in files
        ids, seqs = read_fasta(file)
        df = DataFrame(id=ids, seq=seqs)
        gdf = groupby(df,:seq)
        gdf = combine(gdf,nrow => :count)
        sort!(gdf,:count,rev=true)
        cdf = similar(gdf,1)
        cdf[1,:] = gdf[1,:]
        gdf = gdf[2:end,:]
        clade = 1
        while nrow(gdf) > 0  && clade < max_clades
            dist_cols = PairWise(Hamming,gdf[!,:seq],cdf[!,:seq])
            dist_min = vec(mapslices(minimum,dist_cols,dims=2))
            edf = DataFrame(dm=dist_min)
            gdf = hcat(gdf,edf)
            sort!(gdf,:dm,rev=true)
            max_min_dist = gdf[1,:dm] / length(gdf[1,:seq])
            if verbose println("max_min dist = ", max_min_dist) end
            if max_min_dist < thresh_hold
                println("  breaking ...")
                break
            end
            select!(gdf,Not(:dm))
            push!(cdf,gdf[1,:])
            gdf = gdf[2:end,[:seq,:count]]
            clade += 1
        end
        if verbose println(" $(clade) clades identified in file $(file)") end
        for i in 1:nrow(cdf)
            push!(queries, [replace(basename(file), r".fasta" => "")*"_$(i)"
                            cdf[i,:count]
                            degap(cdf[i,:seq])])
        end
    end
    total_queries = nrow(queries)
    println("BLASTing $(total_queries) query sequences ... ")
    results = DataFrame(sample=String[],score=Float64[],accession=String[],description=String[])
    println(queries[:,:query],queries[:,:sample],max_waits)
    blast_arr = WebBLAST(queries[:,:query],query_names=queries[:,:sample],
                            max_waits=max_waits, num_hits=1)
    if blast_arr != nothing
        rdf = flatten_to_dataframe(blast_arr)
        for (sample, identity, query_len, hit_acc, hit_def) in zip(rdf[:,"Query-def"], rdf[:,:Hsp_identity], rdf[:,"Query-len"], rdf[:,:Hit_accession], rdf[:,:Hit_def] )
            score = round(tryparse(Int64, identity) / tryparse(Int64, query_len),digits=3)
            if verbose
                println("sample = $(sample), score = $score, acc = $(hit_acc)")
                println("$(hit_def)")
            end
            push!(results, [sample score hit_acc hit_def])
        end
    end
    return(results)
end

function format_tbl(df)
    # df = CSV.read(file, DataFrame)
    tbl_html = DataFrames.repr("text/html", df)
    tbl_fmt = replace(tbl_html,
        "<table class=\"data-frame\">" => "<table class=\"table table-striped\";>")
    tbl_fmt = replace(tbl_fmt, r"<p>.+columns</p>" => "") #remove display description
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Int64\">Int64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Int64\">Int64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Float64\">Float64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<tr><th></th></tr>" => "") #remove type rows
    # tbl_fmt = replace(tbl_fmt, "<tr><th></th><th title=\"String\">String</th><th title=\"Int64\">Int64</th></tr>" => "") #remove types
    # tbl_fmt = replace(tbl_fmt, r"<thead>.+</thead>" => "") #remove heading row
    tbl_fmt = replace(tbl_fmt, r"<th></th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th>.</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th>..</th>" => "") #remove row count column
    return tbl_fmt
end


# sample = snakemake.wildcards["sample"]
dataset = snakemake.wildcards["dataset"]

blast_df = DataFrame()
br_tbl = "no blast reports generated, try using SnakeBlast ..."
files = [ "postproc/$(dataset)/$(sample)/$(sample)-blast.csv" for sample in snakemake.params["SAMPLES"] ]
for file in files
    global blast_df = vcat(blast_df, CSV.read(file, DataFrame))
end
CSV.write(snakemake.output[2], blast_df)
    if nrow(blast_df) > 0
        br_tbl = format_tbl(blast_df)
    else
        br_tbl = " BLAST searches returned no hits. "
    end


html_str = """
<html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
        <style>
            body{ margin:0 100; background:white; }
        </style>
        <style type="text/css">
            .container{
            display: flex;
            }
            .col{
                flex: 1;
                text-align: center;
                align-self: center;
            }
            img {
                margin: auto;
            }
        </style>
    </head>
    <body>
        <div style="max-width: 700px;">
            <h1>PorpidPostproc, Blast results</h1>
            <h4><i>Version: $(snakemake.params["VERSION"])</i></h4>
            <h4><i>Commit ID: $(snakemake.params["COMMIT"])</i></h4>
            
            <h3> BLAST results for clades identified in each sample in $(dataset) </h4>
            $(br_tbl)
            
            <h4>blast parameters:</h4>
            <ul>
              <li> thresh_hold = $(snakemake.params["thresh_hold"]) </li>
              <li> max_clades = $(snakemake.params["max_clades"]) </li>
              <li> max_waits = $(snakemake.params["max_waits"]) </li>
            </ul>
        </div>
    </body>
 </html>
""";

open(snakemake.output[1],"w") do io
    write(io,html_str)
end

