#!/bin/bash
./erne-create --fasta toy.reference.fa --output-prefix out
./erne-map --reference out.ebh --query1 toy.long.reads.fastq --output output.long.reads.sam
./erne-map --reference out.ebh --query1 toy.short.reads.fastq --output output.short.reads.sam

