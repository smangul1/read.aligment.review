#!/bin/bash
erne-create --fasta multimapped.reference.fa --output-prefix ref
erne-map --reference ref.ebh --query1 shortReads.fq --output reads.toy.example.erne.bam
