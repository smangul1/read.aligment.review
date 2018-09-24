create index: 

subread-buildindex -o <output dir> reference.fasta

map reads:

subread-align -t 1 -i <index fies> -r reads.toy.example.fasta -o output.bam
