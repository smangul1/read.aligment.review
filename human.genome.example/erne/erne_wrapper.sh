#!/bin/bash
ernecreate=/u/home/v/victorx/changes.read.aligment.review/tools/erne-2.1.1-linux/bin/erne-create
ernemap=/u/home/v/victorx/changes.read.aligment.review/tools/erne-2.1.1-linux/bin/erne-map
read1=/u/home/v/victorx/vic_scratch/SRR043348_1.fastq
read2=/u/home/v/victorx/vic_scratch/SRR043348_2.fastq
genome=/u/home/s/serghei/project/Homo_sapiens/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

#Output is a bam file by default
$ernecreate --fasta $genome --output-prefix /u/home/v/victorx/vic_scratch/ToolsDone/erne/out
$ernemap --reference  /u/home/v/victorx/vic_scratch/ToolsDone/erne/out.ebh --query1 $read1 --query2 $read2 --output /u/home/v/victorx/vic_scratch/ToolsDone/erne/erne.out.bam
cp /u/home/v/victorx/changes.read.aligment.review/tools/erne-2.1.1-linux/bin/erne_wrapper.sh /u/home/v/victorx/vic_scratch/ToolsDone/erne
