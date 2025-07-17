# poretrim.jl

Julia version 1.10.5

`poretrim.jl` is a julia utility for trimming Nanopore `fastq.gz` datasets
using primers known to reside in the dataset.

usage:
```
julia poretrim.jl infile.fastq.gz outfile.fastq.gz fwd_primer rev_primer
```

for example:
```
julia poretrim.jl infile.fastq.gz outfile.fastq.gz TAGGCATCTCCTATGGCAGGAAGAA CCGCTCCGTCCGACGACTCACTATA
```

the tool will read sequence records, search for the forward and reverse primers,
and output a sequence record with the sequence trimmed upto but not including the primers.
If the primer pair is not found the reverse complement sequence is searched and
if still not found the sequence is discarded.

The tool reports every 10K sequences found and at the end produces a histogram 
displaying various sequence counts such as the one below.

![trimcounts](Pool1_Nanopore_trimmed_trimcounts.png)

When you first run the tool, julia will use the Manifest and Porject files to install
the required libraries locally. Thereafter run times should improve.

This is still a work in progress

