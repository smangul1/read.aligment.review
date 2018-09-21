#!/bin/bash
cushaw3=/u/home/v/victorx/changes.read.aligment.review/tools/cushaw3-v3.0.3/cushaw3
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fastq
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fastq
genome=/u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

$cushaw3 index -p genome $genome
$cushaw3 align -r genome -q $read1 $read2 -t 12 | $samtools view -bS - > cushaw3.out.bam
 
