#!/bin/bash

cp short.multimapped.reads.fastq reads.toy.example.fastq
cp short.multimapped.reads.fasta reads.toy.example.fasta
cp long.multimapped.reads.fastq reads.toy.example.long.fastq
cp long.multimapped.reads.fasta reads.toy.example.long.fasta
cp multimapped.reference.short.fasta reference.fasta
cp multimapped.reference.long.fasta reference.long.fasta

segemehl.x -x reference.idx -d reference.fasta
segemehl.x -i reference.idx -d reference.fasta -q reads.toy.example.fasta > results/multimapped.segemehl.sam

smalt index -k 14 -s 8 reference reference.fasta
smalt map -o results/multimapped.smalt.sam reference reads.toy.example.fastq

2bwt-builder reference.fasta
soapsplice -d reference.fasta.index -1 reads.toy.example.fastq -o results/multimapped.soapsplice -f2

splazers reference.fasta reads.toy.example.fasta -o results/multimapped.splazers.result

subread-buildindex -o reference reference.fasta
subread-align -t 1 -i reference -r reads.toy.example.fastq -o results/multimapped.subread.sam --SAMoutput

gmapper-ls reads.toy.example.fasta reference.fasta > results/multimapped.shrimp.sam

nanoblaster -C10 -r reference.long.fasta -i reads.toy.example.long.fastq -o results/multimapped.nanoblaster

micro_razers reference.fasta reads.toy.example.fasta -o results/multimapped.micro_razers.sam

nucmer -p results/multimapped.mummer reference.fasta reads.toy.example.fasta

bowtie-build reference.fasta reference
bowtie -S reference reads.toy.example.fastq > results/multimapped.bowtie.sam    

bowtie2-build reference.fasta reference.fasta 
bowtie2 -x reference.fasta -U reads.toy.example.fastq > results/multimapped.bowtie2.sam

bwa index reference.fasta
bwa mem reference.fasta reads.toy.example.fastq > results/multimapped.bwa.sam

erne-create --fasta reference.fasta --output-prefix ref
erne-map --force-standard --sam --reference ref.ebh --query1 reads.toy.example.fastq --output results/multimapped.erne.sam

hisat2-build reference.fasta reference
#Can use either fasta (-f) or fastq file (-q) for reads
hisat2 -q -x reference -U reads.toy.example.fastq > results/multimapped.hisat2.sam

ngm -r reference.fasta -q reads.toy.example.fastq -o results/multimapped.ngm.sam 

minimap2 -a reference.fasta reads.toy.example.fastq  > results/multimapped.minimap2.sam

cp reads.toy.example.fastq reads.toy.example.fq
cp reference.fasta reference.fa
minialign -xont.1dsq reference.fa reads.toy.example.fq > results/multimapped.minialign.sam

lastdb -uNEAR -R01 reference reference.fasta
lastal -Q1 reference reads.toy.example.fastq | last-split> results/multimapped.last.maf

ngmlr -r reference.fasta -q reads.toy.example.fastq -o results/multimapped.ngmlr.sam

lordfast --index reference.fasta
lordfast --search reference.fasta --seq reads.toy.example.fastq > results/multimapped.lordfast.sam

graphmap align -r reference.fasta -d reads.toy.example.fastq > results/multimapped.graphmap.sam

rm reads.toy.example.fq
rm ref*
rm reads.toy.example.fastq
rm reads.toy.example.fasta
rm reads.toy.example.long.fastq
rm reads.toy.example.long.fasta