ENV["MPLBACKEND"] = "Agg"
using PORPIDpipeline, Base64, CSV, DataFrames, NextGenSeqUtils, FASTX, WebBlast



function format_tbl(df)
    # df = CSV.read(file, DataFrame)
    tbl_html = DataFrames.repr("text/html", df)
    tbl_fmt = replace(tbl_html,
        "<table class=\"data-frame\">" => "<table class=\"table table-striped\";>")
    tbl_fmt = replace(tbl_fmt, r"<p>.+columns</p>" => "") #remove display description
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String15\">String15</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String15\">String31</th>" => "")
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



sample = snakemake.wildcards["sample"]
dataset = snakemake.wildcards["dataset"]

blast_results = DataFrame()
br_tbl = "no blast hits for $(sample)..."

max_clades = snakemake.params["max_clades"]
max_waits = snakemake.params["max_waits"]
thresh_hold = snakemake.params["thresh_hold"]
blast_files = []
# check for panel rejects
ids, seqs = read_fasta(snakemake.input[2])
if (length(seqs)>0)
    blast_files = [snakemake.input[1], snakemake.input[2]]
else
    blast_files = [snakemake.input[1]]
end
t1 = time()
blast_results = get_blast_results(blast_files,thresh_hold=thresh_hold,
                         max_clades=max_clades,max_waits=max_waits)

#    reject_blast_results = DataFrame()
#    reject_blast_results = get_blast_results(snakemake.input[2], thresh_hold=thresh_hold, max_clades=max_clades, max_waits=max_waits)
#    blast_results = vcat(blast_results,reject_blast_results)

if nrow(blast_results) > 0
    br_tbl = format_tbl(blast_results)
else
    br_tbl = " BLAST search returned no hits for $(sample). "
    exit(1)
end
CSV.write(snakemake.output[2], blast_results)
t2 = time()
exec_time = round( (t2-t1), digits=0)
println("$(sample), BLAST search took $(exec_time) seconds.")
br_tbl = br_tbl * "<p> BLAST search for $(sample) took $(exec_time) seconds </p>"

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
            <h2>PORPIDpipeline report</h2>
            <h3><i>Sample: $(sample)</i></h3>
            <h3><i>Dataset: $(dataset)</i></h3>
            <h4><i>Version: $(snakemake.params["VERSION"])</i></h3>
            <h4><i>Commit ID: $(snakemake.params["COMMIT"])</i></h3>
            
            <p>
            <h3> BLAST results for clades identified in $(sample) </h3>
            $(br_tbl)
        </div>
    </body>
 </html>
""";

open(snakemake.output[1],"w") do io
    write(io,html_str)
end
