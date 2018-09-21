#!/bin/bash
./hobbes-index --sref toy.reference.fa -i toy.reference.hix -g 11
./hobbes -q toy.short.reads.fq --sref toy.reference.fa -i toy.reference.hix -a --indel -v 5 --mapout output.short.reads.sam
./hobbes -q toy.long.reads.fq --sref toy.reference.fa -i toy.reference.hix -a --indel -v 5 --mapout output.long.reads.sam
