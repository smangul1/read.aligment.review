#!/bin/bash
erne-create --fasta reference.fasta --output-prefix ref
erne-map --force-standard --reference ref.ebh --query1 reads.toy.example.fastq --output reads.toy.example.erne.bam