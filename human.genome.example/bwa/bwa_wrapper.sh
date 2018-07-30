#!/bin/bash
source /u/local/Modules/default/init/modules.sh
module load bwa
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fastq
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fastq
genome=/u/home/v/victorx/vic_scratch/genome.fa
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools

bwa index $genome 
bwa mem $genome $read1 $read2| $samtools view -bS - > /u/home/v/victorx/vic_scratch/ToolsDone/bwa
cp /u/home/v/victorx/changes.read.aligment.review/tools/bwa/bwa_wrapper.sh /u/home/v/victorx/vic_scratch/ToolsDone/bwa
mv /u/home/v/victorx/changes.read.aligment.review/tools/bwa/annular* /u/home/v/victorx/vic_scratch/ToolsDone/bwa