#!/bin/bash
#Argument 1 is the file name to copy the fastq and fasta files to
cp ./toy.test.data/toy.long_reads.fq $1/toy.long.reads.fq
cp ./toy.test.data/toy.short_reads.fq $1/toy.short.reads.fq
cp ./toy.test.data/toy.reference.fa $1/toy.reference.fa


