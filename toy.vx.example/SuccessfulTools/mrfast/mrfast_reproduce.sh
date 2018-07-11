#!/bin/bash
./mrfast --index toy.reference.fa
./mrfast --search toy.reference.fa --seq toy.short.reads.fq -o output.short.reads.sam
/mrfast --search toy.reference.fa --seq toy.long.reads.fq -o output.long.reads.sam
