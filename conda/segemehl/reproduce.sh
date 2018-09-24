create index:

segemehl.x -x reference.idx -d reference.fasta

map reads:

segemehl.x -i reference.idx -d reference.fasta -q reads.toy.example.fasta > output.map
