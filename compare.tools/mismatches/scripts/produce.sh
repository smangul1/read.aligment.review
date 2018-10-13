#!/bin/bash
#Test.sh takes an argument for the number of mismatches desired
./newtestlong.sh $1
#./newtestshort.sh $1
python mismatch_maker_short.py
python mismatch_maker_long.py
bowtie2 -x reference.fasta -U short.reads.mismatches.fastq > short.reads.mismatches.bowtie2.sam
bowtie2 -x reference.fasta -U long.reads.mismatches.fastq > long.reads.mismatches.bowtie2.sam
cp short.reads.mismatches.fastq ..
cp long.reads.mismatches.fastq ..
cp short.reads.mismatches.bowtie2.sam ..
cp long.reads.mismatches.bowtie2.sam ..