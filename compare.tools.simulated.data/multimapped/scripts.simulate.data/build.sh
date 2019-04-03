#!/bin/bash
#Build the reference and multimapped reads
./createReference.sh
./createShortReads.sh
./createLongReads.sh
python buildMultiRef.short.py
python buildMultiRef.long.py
sed -n '1~4s/^@/>/p;2~4p' shortReads.fq > ../short.multimapped.reads.fasta
sed -n '1~4s/^@/>/p;2~4p' longReads.fq > ../long.multimapped.reads.fasta
mv multimapped.reference.short.fa ../multimapped.reference.short.fasta
mv multimapped.reference.long.fa ../multimapped.reference.long.fasta 
mv shortReads.fq ../short.multimapped.reads.fastq
mv longReads.fq ../long.multimapped.reads.fastq
