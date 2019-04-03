#!/bin/bash

segemehl.x -x reference.idx -d reference.fasta
segemehl.x -i reference.idx -d reference.fasta -q reads.toy.example.fasta > results/short.unmapped.segemehl.sam

smalt index -k 14 -s 8 reference reference.fasta
smalt map -o results/short.unmapped.smalt.sam reference reads.toy.example.fastq

2bwt-builder reference.fasta
soapsplice -d reference.fasta.index -1 reads.toy.example.fasta -o results/short.unmapped.soapsplice.sam -f2

splazers reference.fasta reads.toy.example.fasta -o results/short.unmapped.splazers.result

subread-buildindex -o reference reference.fasta
subread-align -t 1 -i reference -r reads.toy.example.fastq -o results/short.unmapped.subread.sam --SAMoutput


gmapper-ls reads.toy.example.fasta reference.fasta > results/short.unmapped.shrimp.sam

nanoblaster -C10 -r reference.fasta -i reads.toy.example.long.fastq -o results/short.unmapped.nanoblaster

micro_razers reference.fasta reads.toy.example.short.fasta -o results/short.unmapped.micro_razers.sam

nucmer -p results/short.unmapped.mummer reference.fasta reads.toy.example.fasta

bowtie-build reference.fasta reference
bowtie -S reference reads.toy.example.short.fastq > results/short.unmapped.bowtie.sam    

bowtie2-build reference.fasta reference.fasta 
bowtie2 -x reference.fasta -U reads.toy.example.fastq > results/short.unmapped.bowtie2.sam

bwa index reference.fasta
bwa mem reference.fasta reads.toy.example.fastq > results/short.unmapped.bwa.sam

erne-create --fasta reference.fasta --output-prefix ref
erne-map --force-standard --sam --reference ref.ebh --query1 reads.toy.example.short.fastq --output results/short.unmapped.erne.sam

hisat2-build reference.fasta reference
#Can use either fasta (-f) or fastq file (-q) for reads
hisat2 -q -x reference -U reads.toy.example.fastq > results/short.unmapped.hisat2.sam

ngm -r reference.fasta -q reads.toy.example.fastq -o results/short.unmapped.ngm.sam 

minimap2 -a reference.fasta reads.toy.example.fastq  > results/short.unmapped.minimap2.sam

cp reads.toy.example.fastq reads.toy.example.fq
cp reference.fasta reference.fa
minialign -xont.1dsq reference.fa reads.toy.example.fq > results/short.unmapped.minialign.sam

lastdb -uNEAR -R01 reference reference.fasta
lastal -Q1 reference reads.toy.example.fastq | last-split> results/short.unmapped.last.maf

ngmlr -r reference.fasta -q reads.toy.example.fastq -o results/short.unmapped.ngmlr.sam

lordfast --index reference.fasta
lordfast --search reference.fasta --seq reads.toy.example.fastq > results/short.unmapped.lordfast.sam

graphmap align -r reference.fasta -d reads.toy.example.fastq > results/short.unmapped.graphmap.sam

rm reads.toy.example.fq
cp reference.fasta copy.fa
rm ref*
mv copy.fa reference.fasta