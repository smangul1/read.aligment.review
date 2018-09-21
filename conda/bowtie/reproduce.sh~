#!/bin/bash
bowtie-build reference.fasta reference
bowtie -S reference reads.toy.example.fastq | samtools view -bS - > reads.toy.example.bowtie.bam  

#reads limited to 1024 characters