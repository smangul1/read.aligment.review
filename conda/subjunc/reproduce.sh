create index: 

subread-buildindex -o <output folder> reference.fasta

map reads:

subjunc -T 5 -i <index file> -r reads.toy.example.fasta -o output.bam

