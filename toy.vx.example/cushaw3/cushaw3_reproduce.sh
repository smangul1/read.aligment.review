#!/bin/bash

#Constructs the BWT and FM-index for reference genome toy.reference.fa
cushaw3 index -p toy.reference toy.reference.fa
#Aligner
#-r bwt_file_base -> reference genome that was indexed in previous step
#-f set of reads to be aligned in fa/fq format 
#-t numThreads 
#-o sam_file
./cushaw3 align -r toy.reference -f toy.short.reads.fq -t 12 -o output.short.reads.sam
./cushaw3 align -r toy.reference -f toy.long.reads.fq -t 12 -o output.long.reads.sam




