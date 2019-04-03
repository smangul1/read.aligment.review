#!/bin/bash

cp shortest.read.length.fastq reads.toy.example.fastq
cp shortest.read.length.fasta reads.toy.example.fasta 

blasr --header reference.fasta reads.toy.example.fasta > results/shortest.blasr.result

snap-aligner index reference.fasta index-dir
snap-aligner single index-dir reads.toy.example.fastq -o results/shortest.snap.sam

rmap reads.toy.example.fastq -c reference.fasta -o results/shortest.rmap.sam

segemehl.x -x reference.idx -d reference.fasta
segemehl.x -i reference.idx -d reference.fasta -q reads.toy.example.fasta > results/shortest.segemehl.sam

smalt index -k 14 -s 8 reference reference.fasta
smalt map -o results/shortest.smalt.sam reference reads.toy.example.fastq

#had to manually test, it cannot align a read 10 bps long
2bwt-builder reference.fasta
soapsplice -d reference.fasta.index -1 soapsplice.reads.fastq -o results/shortest.soapsplice -f2

splazers reference.fasta reads.toy.example.fasta -o results/shortest.splazers.result

subread-buildindex -o reference reference.fasta
subread-align -t 1 -i reference -r reads.toy.example.fastq -o results/shortest.subread.sam --SAMoutput

gmapper-ls reads.toy.example.fasta reference.fasta > results/shortest.shrimp.sam

nanoblaster -C10 -r reference.fasta -i reads.toy.example.long.fastq -o results/shortest.nanoblaster

micro_razers reference.fasta reads.toy.example.fasta -o results/shortest.micro_razers.sam

nucmer -p results/shortest.mummer reference.fasta reads.toy.example.fasta

bowtie-build reference.fasta reference
bowtie -S reference reads.toy.example.fastq > results/shortest.bowtie.sam    

bowtie2-build reference.fasta reference.fasta 
bowtie2 -x reference.fasta -U reads.toy.example.fastq > results/shortest.bowtie2.sam

bwa index reference.fasta
bwa mem reference.fasta reads.toy.example.fastq > results/shortest.bwa.sam

erne-create --fasta reference.fasta --output-prefix ref
erne-map --force-standard --sam --reference ref.ebh --query1 reads.toy.example.fastq --output results/shortest.erne.sam

hisat2-build reference.fasta reference
#Can use either fasta (-f) or fastq file (-q) for reads
hisat2 -q -x reference -U reads.toy.example.fastq > results/shortest.hisat2.sam

ngm -r reference.fasta -q reads.toy.example.fastq -o results/shortest.ngm.sam 

minimap2 -a reference.fasta reads.toy.example.fastq  > results/shortest.minimap2.sam

cp reads.toy.example.fastq reads.toy.example.fq
cp reference.fasta reference.fa
minialign -xont.1dsq reference.fa reads.toy.example.fq > results/shortest.minialign.sam

lastdb -uNEAR -R01 reference reference.fasta
lastal -Q1 reference reads.toy.example.fastq | last-split> results/shortest.last.maf

ngmlr -r reference.fasta -q reads.toy.example.fastq -o results/shortest.ngmlr.sam

lordfast --index reference.fasta
lordfast --search reference.fasta --seq reads.toy.example.fastq > results/shortest.lordfast.sam

graphmap align -r reference.fasta -d reads.toy.example.fastq > results/shortest.graphmap.sam

rm reads.toy.example.fq
cp reference.fasta copy.fa
rm ref*
mv copy.fa reference.fasta
rm reads.toy.example.fastq
rm reads.toy.example.fasta