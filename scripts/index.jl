ENV["MPLBACKEND"] = "Agg"
using Base64, CSV, DataFrames, NextGenSeqUtils, FASTX, WebBlast

function get_image_str(file)
    fig_str = open(file,"r") do io
        read(io)
    end
    return Base64.base64encode(fig_str)
end



function format_tbl(df)
    # df = CSV.read(file, DataFrame)
    tbl_html = DataFrames.repr("text/html", df)
    tbl_fmt = replace(tbl_html,
        "<table class=\"data-frame\">" => "<table class=\"table table-striped\";>")
    tbl_fmt = replace(tbl_fmt, r"<p>.+columns</p>" => "") #remove display description
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String\">String</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String31\">String31</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"String15\">String15</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Any\">Any</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Int64\">Int64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Int64\">Int64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Float64\">Float64</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Bool\">Bool</th>" => "")
    tbl_fmt = replace(tbl_fmt, r"<th title=\"Bool\">Bool</th>" => "")
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

# qc_fig = get_image_str(snakemake.input[1]);
# qc_tbl = format_tbl(snakemake.input[2]);
# mds_fig = get_image_str(snakemake.input[3]);
# phylo_fig = get_image_str(snakemake.input[4]);

qual_df = CSV.read("porpid/$(dataset)/quality_report.csv", DataFrame)
qual_tbl = format_tbl(qual_df)

demux_df = CSV.read("porpid/$(dataset)/demux_report.csv", DataFrame)
demux_tbl = format_tbl(demux_df)

####### make a table of ALL read counts
seq_counts_df = DataFrame(Sample = [], Porpid_Seqs = [], Rejected_Min_Ag = [], Rejected_Panel = [], Rejected_Seqs = [], Final_Seqs = [])
for sample in sort(snakemake.params["SAMPLES"])
    p_seqs, p_seq_names = read_fasta("porpid/$(dataset)/consensus/$(sample).fasta")    #porpid sequences
    r_seqs, r_seq_names = read_fasta("postproc/$(dataset)/$(sample)/$(sample).fasta.rejected.fasta") #rejected sequences
    f_seqs, f_seq_names = read_fasta("postproc/$(dataset)/$(sample)/$(sample).fasta") #final sequences
    sample_reject_df = CSV.read("postproc/$(dataset)/$(sample)/$(sample).fasta.rejected.csv", DataFrame) #reject split
    r_ma_seq_number = sample_reject_df[1,"count"]
    r_pan_seq_number = sample_reject_df[2,"count"]
    p_seq_number = length(p_seqs)
    r_seq_number = length(r_seqs)
    f_seq_number = length(f_seqs)
    push!(seq_counts_df, [sample, p_seq_number, r_ma_seq_number, r_pan_seq_number, r_seq_number, f_seq_number])
end
seq_counts_df[!, :Porpid_Seqs] = convert.(Int, seq_counts_df[:, :Porpid_Seqs])
seq_counts_df[!, :Rejected_Min_Ag] = convert.(Int, seq_counts_df[:, :Rejected_Min_Ag])
seq_counts_df[!, :Rejected_Panel] = convert.(Int, seq_counts_df[:, :Rejected_Panel])
seq_counts_df[!, :Rejected_Seqs] = convert.(Int, seq_counts_df[:, :Rejected_Seqs])
seq_counts_df[!, :Final_Seqs] = convert.(Int, seq_counts_df[:, :Final_Seqs])

#create final table with sequence number and reads per template for porpid seqs
joined_df = innerjoin(seq_counts_df, demux_df, on = :Sample)
joined_df = rename!(joined_df,:Count => :Read_Count) #change Counts column name to Read_Count
joined_df[!, :Reads_per_Porpid_Seq] = joined_df[!, :Read_Count] ./ joined_df[!, :Porpid_Seqs]
joined_df = select(joined_df, [:Sample, :Reads_per_Porpid_Seq], :Porpid_Seqs, :Rejected_Min_Ag, :Rejected_Panel, :Rejected_Seqs, :Final_Seqs)
joined_df_tbl = format_tbl(joined_df)
CSV.write(snakemake.output[3], joined_df)


contam_df = CSV.read("porpid/$(dataset)/contam_report.csv", DataFrame)
# contam_df = contam_df[contam_df[:,:discarded],:]
contam_df = filter(row -> row.discarded == "true", contam_df)
contam_df = select(contam_df, :sample, :sequence_name, :nearest_nonself_variant, :nearest_nonself_distance)
contam_tbl = format_tbl(contam_df)
contam_tbl = replace(contam_tbl,
     "<table class=\"table table-striped\";>" => "<table class=\"table table-striped\";  style=\"font-size: 8px\";>")
num_discarded = nrow(contam_df)
if num_discarded == 0
    contam_tbl = ""
end

max_suspects = 10
contam_suspect_df = CSV.read("porpid/$(dataset)/contam_suspect.csv", DataFrame)
num_suspected = nrow(contam_suspect_df)
if nrow(contam_suspect_df) > max_suspects
    contam_suspect_df = contam_suspect_df[1:max_suspects,:]
end
contam_suspect_tbl = format_tbl(contam_suspect_df)
contam_suspect_tbl = replace(contam_suspect_tbl,
     "<table class=\"table table-striped\";>" => "<table class=\"table table-striped\";  style=\"font-size: 8px\";>")
if num_suspected == 0
    contam_suspect_tbl = ""
end

# contam_ptr = "click <a href=$(dataset)-contam.html> here </a> for contam report."


reject_df = CSV.read("porpid/$(dataset)/reject_report.csv", DataFrame)
reject_tbl = format_tbl(reject_df);

demux_dict = Dict(Pair.(demux_df.Sample, demux_df.Count))


br_tbl = "click <a href=$(dataset)-blast.html> here </a> for blast reports."

# get parameters for parameter reporting



html_str_hdr = """
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
"""
html_str = html_str_hdr * """
    <body>
        <div style="max-width: 700px;">
            <h1>PORPIDpipeline</h1>
            <h4><i>Version: $(snakemake.params["VERSION"])</i></h4>
            <h4><i>Commit ID: $(snakemake.params["COMMIT"])</i></h4>
            
            <h1>Dataset: $(dataset)</h1>
            
            <h3>demux quality report:<h3>
            $(qual_tbl)
            
            <h3>demux rejection report:</h3>
            $(reject_tbl)
            
            <h4>demux parameters:</h4>
            <ul>
              <li> chunk_size = $(snakemake.params["chunk_size"]) </li>
              <li> error_rate = $(snakemake.params["error_rate"]) </li>
              <li> min_length = $(snakemake.params["min_length"]) </li>
              <li> max_length = $(snakemake.params["max_length"]) </li>
            </ul>
"""
                     
if num_discarded == 0
    contam_hdr = "$(num_discarded) sequences are classified as contaminates."
elseif num_discarded == 1
    contam_hdr = "$(num_discarded) sequence is classified as contaminated and has been discarded."
else
    contam_hdr = "$(num_discarded) sequences are classified as contaminated and have been discarded."
end

suspected = nrow(contam_suspect_df)
if num_suspected == 1
    contam_suspect_hdr = "$(num_suspected) sequence has been flagged as a possible contaminate."
else
    contam_suspect_hdr = "$(num_suspected) sequences have been flagged as possible contaminates, here are the top $(suspected)."
end
            
html_str = html_str * """
                <h3>contamination report:</h3>
                <h4> contam filter toggle switch is
                $(snakemake.params["contam_toggle"])</h4>
                <h4>discarded contaminated sequences:</h4>
                $(contam_hdr) <br>
                $(contam_tbl) <br>
                A full report is available in <tt>contam-report.csv</tt>
                <h4>suspected contaminates:</h4>
                $(contam_suspect_hdr) <br>
                $(contam_suspect_tbl)
                A full list of suspected contaminates is available in <tt>contam-suspect.csv</tt>
                <h4>contam parameters:</h4>
                <ul>
                  <li> contam_toggle = $(snakemake.params["contam_toggle"]) </li>
                  <li> cluster_thresh = $(snakemake.params["cluster_thresh"]) </li>
                  <li> error_rate = $(snakemake.params["error_rate"]) </li>
                  <li> proportion_thresh = $(snakemake.params["proportion_thresh"]) </li>
                  <li> dist_thresh = $(snakemake.params["dist_thresh"]) </li>
                </ul>
""";

html_str = html_str * """
            <h3>sample reports:</h3>
            <div class=\"data-frame\"><table class=\"table table-striped\";><thead>
            <tr>
              <th>Sample</th>
              <th style=\"text-align:center\">Count</th>
              <th style=\"text-align:left\">LR%</th>
            </tr></thead>
""";

for sample in sort(snakemake.params["SAMPLES"])
    qc_bins = CSV.read("postproc/$(dataset)/$(sample)/$(sample)_qc_bins.csv",DataFrame)
    success = 0
    nr = nrow(qc_bins)
    nc = ncol(qc_bins)
    if qc_bins[nr,1] == "likely_real"
        success = floor(Int, 100 * qc_bins[nr,nc] / sum(qc_bins[:,nc]))
    end
    
    global html_str = html_str * "<tr> <td><a href=$(sample)/$(sample)-report.html target=blank> $(sample) </a> </td> <td style=\"text-align:center\"> $(demux_dict[sample]) </td>  <td> $(success)% </td></tr>"
end

html_str = html_str * """
    </table></div>
    <h4>parameters:</h4>
    <ul>
      <li> fs_thresh = $(snakemake.params["fs_thresh"]) </li>
      <li> lda_thresh = $(snakemake.params["lda_thresh"]) </li>
      <li> agreement_thresh = $(snakemake.params["agreement_thresh"]) </li>
      <li> panel_thresh = $(snakemake.params["panel_thresh"]) </li>
    </ul>
    <h3>sequence counts:</h3>
    $(joined_df_tbl)
    Summary of sequence output from porpid, those that were rejected for being too disimilar to
    the panel file, and the final sequence count after filtering. Reads per porpid sequence
    can be used to compare average depth across different samples. 
""";



html_str = html_str * """
            <h4> BLAST results for clades identified in each sample </h4>
            $(br_tbl)
        </div>
    </body>
 </html>
""";

open(snakemake.output[1],"w") do io
    write(io,html_str)
end

html_str = html_str_hdr * """
<body>
<h3> Blast results ... </h3>

No blast results available yet, to get blast results run

<code>
snakemake -s SnakeBlast -k --rerun-incomplete -j1
</code>
</body>
"""

open(snakemake.output[2],"w") do io
    write(io,html_str)
end





