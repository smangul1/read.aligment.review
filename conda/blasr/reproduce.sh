#!/bin/bash
blasr --header reference.fasta reads.toy.example.fasta > reads.toy.example.blasr

#blasr has a sam and bam option: blasr reads.fasta ref.fasta --bam --out out.bam   However, the command gives ERROR, can not convert non-pacbio reads to pbbam record'?   According to the Github FAQ, the fasta input must follow a specific convention
