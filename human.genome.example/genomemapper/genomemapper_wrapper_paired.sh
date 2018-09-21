#!/bin/bash
gmindex=/u/home/v/victorx/vic_scratch/genomemapper-0.4.4/gmindex
genomemapper=/u/home/v/victorx/vic_scratch/genomemapper-0.4.4/genomemapper
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fastq
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fastq
genome=/u/home/v/victorx/vic_scratch/genome.fa
#Parameters were taken from the quick start of the manual
#-s Hash size/ seed size of 11
#-M 4 mismatches -G 2 Gaps -E Total of 4 edit operations
$gmindex -i $genome -s 11
$genomemapper -i $genome -q $read1 $read2 -o /u/home/v/victorx/vic_scratch/ToolsDone/genomemapper2/genomemapper.pairedread.output.list

cp /u/home/v/victorx/vic_scratch/genomemapper-0.4.4/genomemapper_wrapper2.sh /u/home/v/victorx/vic_scratch/ToolsDone/genomemapper2/genomemapper_wrapper2.sh

#gmindex will create genome.cid, genome.mfd, genome.mrc genome.ta