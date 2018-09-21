#!/bin/bash
./MosaikBuild -fr toy.reference.fa -oa toy.reference.dat
./MosaikBuild -q toy.short.reads.fq -st illumina -out toy.short.reads.mkb
cp ../src/networkFile/2.1.78.pe.ann ./pe.ann
cp ../src/networkFile/2.1.78.pe.ann ./se.ann
./MosaikAligner -in toy.short.reads.mkb -out output.short.reads.mka -ia toy.reference.dat -annpe pe.ann -annse se.ann
