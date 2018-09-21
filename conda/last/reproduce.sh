#!/bin/bash
lastdb -uNEAR -R01 reference reference.fasta
lastal -Q1 reference reads.toy.example.fastq | last-split> reads.toy.example.last.maf
