#!/bin/bash

mummer=/u/home/v/victorx/changes.read.aligment.review/tools/mummer-4.0.0beta2/mummer
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fasta
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fasta
genome=/u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

#mummer efficiently locates maximal unique matches between two sequences
#nucmer for the alignment of multiple closely related nucleotide sequences
#nucmer outputs to a delta file which can be read using mummer's view tools
$mummer $genome $read1 $read2 > /u/home/v/victorx/vic_scratch/ToolsDone/mummer/mummer.output

cp /u/home/v/victorx/changes.read.aligment.review/tools/mummer-4.0.0beta2/mummer_wrapper.sh /u/home/v/victorx/vic_scratch/ToolsDone/mummer/mummer_wrapper.sh