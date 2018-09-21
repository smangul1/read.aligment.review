#!/bin/bash
bowtie-build reference.fasta reference
mapsplice.py -c . -x reference -1 reads.toy.example.fastq --bam

#-c is the directory the fasta file is located in
