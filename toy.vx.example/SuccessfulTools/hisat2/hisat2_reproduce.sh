#!/bin/bash
./hisat2-build toy.reference.fa toy.reference
#Can use either fasta (-f) or fastq file (-q) for reads
./hisat2 -q -x toy.reference -U toy.short.reads.fq -S output.short.reads.sam
./hisat2 -q -x toy.reference -U toy.long.reads.fq -S output.long.reads.sam

