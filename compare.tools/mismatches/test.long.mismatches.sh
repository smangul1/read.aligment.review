#!/bin/bash
cp short.reads.mismatches.fastq reads.toy.example.fastq
#skipped fasta reads for now -> micro_razers, mummer
nanoblaster -C10 -r reference.fasta -i reads.toy.example.fastq -o results/short.mismatches.nanoblaster

bowtie-build reference.fasta reference
bowtie -S reference reads.toy.example.fastq > results/short.mismatches.bowtie.sam    

bowtie2-build reference.fasta reference.fasta 
bowtie2 -x reference.fasta -U reads.toy.example.fastq > results/short.mismatches.bowtie2.sam

bwa index reference.fasta
bwa mem reference.fasta reads.toy.example.fastq > results/short.mismatches.bwa.sam

erne-create --fasta reference.fasta --output-prefix ref
erne-map --force-standard --sam --reference ref.ebh --query1 reads.toy.example.fastq --output results/short.mismatches.erne.sam

hisat2-build reference.fasta reference
#Can use either fasta (-f) or fastq file (-q) for reads
hisat2 -q -x reference -U reads.toy.example.fastq > results/short.mismatches.hisat2.sam

ngm -r reference.fasta -q reads.toy.example.fastq -o results/short.mismatches.ngm.sam 

minimap2 -a reference.fasta reads.toy.example.fastq  > results/short.mismatches.minimap2.sam

cp reads.toy.example.fastq reads.toy.example.fq
cp reference.fasta reference.fa
minialign -xont.1dsq reference.fa reads.toy.example.fq > results/short.mismatches.minialign.sam

lastdb -uNEAR -R01 reference reference.fasta
lastal -Q1 reference reads.toy.example.fastq | last-split> results/short.mismatches.last.maf


ngmlr -r reference.fasta -q reads.toy.example.fastq -o results/short.mismatches.ngmlr.sam

lordfast --index reference.fasta
lordfast --search reference.fasta --seq reads.toy.example.fastq > results/short.mismatches.lordfast.sam

graphmap align -r reference.fasta -d reads.toy.example.fastq > results/short.mismatches.graphmap.sam

rm reads.toy.example.fastq
rm reads.toy.example.fq
cp reference.fasta copy.fa
rm reference*
mv copy.fa reference.fasta