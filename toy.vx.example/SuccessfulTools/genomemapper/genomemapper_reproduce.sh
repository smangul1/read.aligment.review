#!/bin/bash
./gmindex -i toy.reference.fa -s 11
./genomemapper -i toy.reference.fa -q toy.long.reads.fq -M 4 -G 2 -E 4 -o output.long.reads.shore
./genomemapper -i toy.reference.fa -q toy.short.reads.fq -M 4 -G 2 -E 4 -o output.short.reads.shore

