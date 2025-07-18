using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using CSV, DataFrames, Plots
using BioSequences, FASTX, YAML, StatsBase
using CodecZlib: GzipDecompressorStream
using CodecZlib: GzipCompressorStream

################################## some functions #########################

function string_reverse_complement(dna_string::String)
    return String(reverse_complement(LongDNA{4}(dna_string)))
end

function primer_peek(seqs::Vector{Any}; L=20, N=30, keep=20 )
    starts=[ s[1:L] for s in seqs ]
    ends=[ s[end-L+1:end] for s in seqs ]
    pot_primers = countmap( vcat(starts,string_reverse_complement.(ends)))
    topN=reverse(sort([(pot_primers[k],k) for k in keys(pot_primers)]))[1:N]
    for i in enumerate(topN)
        println(i)
    end
    # return my_reverse_complement.(sort(my_reverse_complement.([i[2] for i in topN[1:keep]])))
    return topN[1:keep]
end

function findbestlast(query,tol,window)
    best=nothing
    next=findfirst(query,tol,window)
    while !isnothing(next)
        best=next
        next=findnext(query,tol,window,last(next)+1)
    end
    return best
end

function nanopore_trim(in_file,out_file,top,bot;
                            window=150, tol=2,
                            primer_discard_file=nothing,
                            length_discard_file=nothing)
    reader = FASTQ.Reader(GzipDecompressorStream(open(in_file)))
    count=0
    discard_count=0
    top_chop_count=0
    bot_chop_count=0
    both_chop_count=0
    rev_count=0
    short_long_count=0
    seqs=[]
    tps=[]
    bps=[]
    writer = FASTQ.Writer(GzipCompressorStream(open(out_file,"w")))
    if !isnothing(primer_discard_file)
        primer_discard_writer = FASTQ.Writer(GzipCompressorStream(open(primer_discard_file,"w")))
    end
    if !isnothing(length_discard_file)
        length_discard_writer = FASTQ.Writer(GzipCompressorStream(open(length_discard_file,"w")))
    end
    query_top=ApproximateSearchQuery(LongDNA{4}(top))
    query_bot=ApproximateSearchQuery(LongDNA{4}(bot))
    while !eof(reader)
        record = FASTQ.Record()
        read!(reader, record)
        seq=FASTQ.sequence(LongDNA{4},record)
        qual=FASTQ.quality(record)
        count=count+1
        rev=false
        if length(seq) > 2800 && length(seq) < 3800
            top_window=seq[1:window]
            top_chop=findbestlast(query_top, tol, top_window )
            
            bot_window=seq[end-window+1:end]
            bot_chop=findbestlast(query_bot, tol, reverse_complement(bot_window))
            
            if isnothing(top_chop) || isnothing(bot_chop)
                rev_top_window=reverse_complement(seq)[1:window]
                top_chop=findbestlast(query_top, tol, rev_top_window )
                rev_bot_window=reverse_complement(seq)[end-window+1:end]
                bot_chop=findbestlast(query_bot, tol, reverse_complement(rev_bot_window))
                rev=true
            end
            
            if  !isnothing(top_chop)
                top_chop_count = top_chop_count + 1
            end
            if  !isnothing(bot_chop)
                bot_chop_count = bot_chop_count + 1
            end
            if  !isnothing(top_chop) && !isnothing(bot_chop)
                both_chop_count = both_chop_count + 1
                if rev
                    seq=reverse_complement(seq)
                    qual=reverse(qual)
                    rev_count = rev_count + 1
                end
                nam="seq$(both_chop_count)"
                start=first(top_chop)
                start<1 ? start=1 : nothing
                stop=first(bot_chop)
                stop<1 ? stop=1 : nothing
                seq=seq[start:end-stop+1]
                qual=qual[start:end-stop+1]
                # record_fastq=FASTQRecord(nam, seq, qual)
                record_fastq=FASTQRecord(nam, reverse_complement(seq),
                                            reverse(qual))
                write(writer, record_fastq)
                push!(seqs,String(seq))
                push!(tps,start)
                push!(bps,stop)
                if both_chop_count % 10000 == 0
                    @show both_chop_count, count, top_chop_count, bot_chop_count, rev_count, discard_count, short_long_count
                end
            else
                discard_count += 1
                if !isnothing(primer_discard_file)
                    nam="primer_discard_$(discard_count)"
                    record_fastq = FASTQRecord(nam, reverse_complement(seq),reverse(qual))
                    write(primer_discard_writer, record_fastq)
                end
            end
        else
            short_long_count += 1
            if !isnothing(length_discard_file)
                nam="length_discard_$(short_long_count)"
                record_fastq = FASTQRecord(nam, reverse_complement(seq),reverse(qual))
                write(length_discard_writer, record_fastq)
            end
        end
    end
    @show count, top_chop_count, bot_chop_count, both_chop_count, rev_count, short_long_count, discard_count
    close(reader)
    close(writer)
    if !isnothing(primer_discard_file)
        close(primer_discard_writer)
    end
    if !isnothing(length_discard_file)
        close(length_discard_writer)
    end
    return seqs,tps,bps,count,discard_count,short_long_count
end

##################### main script starts here ##############################

if length(ARGS) < 4
    println("************************* usage error **************************")
    println("usage: julia poretrim.jl in_file out_file fwd_primer rev_primer primer_discard_file length_discard_file")
    println("   eg: julia poretrim.jl ../raw-reads/Pool1_nanopore.fastq.gz ./Pool1_nanopore_trimmed.fastq.gz TAGGCATCTCCTATGGCAGGAAGAA CCGCTCCGTCCGACGACTCACTATA primer_discards.fastq.gz length_discards.fastq.gz")
    println("************************* usage error **************************")
    exit()
end

in_file=ARGS[1]
out_file=ARGS[2]
fwd_primer=ARGS[3]
rev_primer=ARGS[4]
length(ARGS) > 4 ? primer_discard_file=ARGS[5] : primer_discard_file=nothing
length(ARGS) > 5 ? length_discard_file=ARGS[6] : length_discard_file=nothing

t1 = time()

window=1000
tol=2
trim_seqs,tps,bps,count,discard_count, short_long_count = nanopore_trim(in_file,out_file,fwd_primer,rev_primer,
                            window=window,tol=tol,primer_discard_file=primer_discard_file,length_discard_file=length_discard_file)

res=primer_peek(trim_seqs,N=6,keep=4,L=min(length(fwd_primer),length(rev_primer)))

t2 = time()

dt=Int( floor( (t2-t1)/60 ) )

histogram(length.(trim_seqs[1:100:end]), label="seq_length",
    xlabel="$(length(trim_seqs)) seqs from $(count) \n with $(short_long_count) size filtered and $(discard_count) primer filtered\n " *
        " with $(window) search window and $(tol) match tolerance in $(dt) minutes",
    title="poretrim filtering \n" * "$(res[1]) \n $(res[2]) \n" )
            
histogram!(tps[1:100:end],label=" 0 + fwd_trim")
histogram!((2*window) .- bps[1:100:end],label="2k - rev_trim")
            
savefig("$(split(basename(out_file),'.')[1])_trimcounts.png")


println("Nanopore trimming took $(t2-t1) seconds.")
