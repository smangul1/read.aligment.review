#!/bin/bash
minimap2 -a reference.fasta reads.toy.example.fastq | samtools view -bS - > reads.toy.example.minimap2.bam
