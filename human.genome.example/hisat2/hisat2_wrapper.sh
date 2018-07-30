#!/bin/bash

hisat2build=/u/home/v/victorx/changes.read.aligment.review/tools/hisat2-2.1.0/hisat2-build
hisat2=/u/home/v/victorx/changes.read.aligment.review/tools/hisat2-2.1.0/hisat2
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fastq
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fastq
genome=/u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

$hisat2build $genome genome
#Can use either fasta (-f) or fastq file (-q) for reads
$hisat2 -q -x genome -1 $read1 -2 $read2 | $samtools view -bS - > hisat2.out.bam
mv hisat2.out.bam /u/home/v/victorx/vic_scratch/hisat2
mv genome* /u/home/v/victorx/vic_scratch/hisat2




