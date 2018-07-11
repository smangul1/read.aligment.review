#!/bin/bash
./segemehl.x -x toy.reference.idx -d toy.reference.fa
#The output doesn't seem to be in a SAM format. 
./segemehl.x -i toy.reference.idx -d toy.reference.fa -q toy.short.reads.fa > output.short.reads.map
./segemehl.x -i toy.reference.idx -d toy.reference.fa -q toy.long.reads.fa > output.long.reads.map