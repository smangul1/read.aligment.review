#!/bin/bash
cp short.read.fastq short.reads.mismatches.fastq

for i in $(seq 1 $1)
do  
    cat short.read.fastq >> short.reads.mismatches.fastq
done
wc -l short.reads.mismatches.fastq