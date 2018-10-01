#!/bin/bash
cp long.read.fastq long.reads.mismatches.fastq

for i in $(seq 1 $1)
do  
    cat long.read.fastq >> long.reads.mismatches.fastq
done
wc -l long.reads.mismatches.fastq