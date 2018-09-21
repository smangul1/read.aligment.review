#!/bin/bash
minialign -xont.1dsq reference.fa reads.toy.example.fq | samtools view -bS - > reads.toy.example.minialign.bam