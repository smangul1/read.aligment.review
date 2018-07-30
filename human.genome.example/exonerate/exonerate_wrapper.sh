#!/bin/bash
exonerate=/u/home/v/victorx/changes.read.aligment.review/tools/exonerate-2.2.0-x86_64/bin/exonerate
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fasta
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fasta
genome=/u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

$exonerate --query $read1 $read2 --target $genome > /u/home/v/victorx/vic_scratch/ToolsDone/exonerate/exonerate.out
mv annular* /u/home/v/victorx/vic_scratch/ToolsDone/exonerate