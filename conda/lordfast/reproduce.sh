#!/bin/bash
lordfast --index reference.fasta
lordfast --search reference.fasta --seq reads.toy.example.fastq | samtools view -bS - > reads.toy.example.lordfast.bam
