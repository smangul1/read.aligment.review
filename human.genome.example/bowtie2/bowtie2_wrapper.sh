#!/bin/bash
source /u/local/Modules/default/init/modules.sh
module load bowtie2
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fastq
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fastq
genome=/u/home/v/victorx/vic_scratch/genome.fa
samtools=/u/home/s/serghei/collab/code/rop.old/tools/samtools

bowtie2-build $genome /u/home/v/victorx/vic_scratch/bowtie2_index
bowtie2 -x /u/home/v/victorx/vic_scratch/bowtie2_index -1 $read1 -2 $read2| $samtools view -bS - > /u/home/v/victorx/vic_scratch/ToolsDone/bowtie2

cp /u/home/v/victorx/changes.read.aligment.review/tools/bowtie2/bowtie2_wrapper.sh /u/home/v/victorx/vic_scratch/ToolsDone/bowtie2

mv /u/home/v/victorx/changes.read.aligment.review/tools/bowtie2/annular* /u/home/v/victorx/vic_scratch/ToolsDone/bowtie2
