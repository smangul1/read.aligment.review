#!/bin/bash
hisat2-build reference.fasta reference
#Can use either fasta (-f) or fastq file (-q) for reads
hisat2 -q -x reference -U reads.toy.example.fastq | samtools view -bS - > reads.toy.example.hisat2.bam
