#!/bin/bash
./mrsfast --index toy.reference.fa
./mrsfast --search toy.reference.fa --seq toy.short.reads.fq -o output.short.reads.sam
/mrsfast --search toy.reference.fa --seq toy.long.reads.fq -o output.long.reads.sam
