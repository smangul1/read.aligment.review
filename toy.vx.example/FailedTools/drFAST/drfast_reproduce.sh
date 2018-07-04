#!/bin/bash
./drfast --index toy.reference.fa
./drfast --search toy.reference.fa --seq toy.short.reads.fq -o output.short.reads.sam
./drfast --search toy.reference.fa --seq toy.long.reads.fq -o output.long.reads.sam
