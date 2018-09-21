#!/bin/bash
drfast --index ref.fasta
drfast --search ref.fasta --seq reads.fastq -o output.sam
