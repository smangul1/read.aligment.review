#!/bin/bash
#createReference.sh, createShortReads.sh  must be run first to get the necessary files used for this script
rm multi.ref.fa
rm multi.read.fq

#Add header to multimapped reference file
head -1 reference.fa >> multi.ref.fa

#Take the first 4 lines of refcopy and duplicate it. Refcopy is reference.fa separated by new lines and without a header
for ((i=1;i<=5;i++));
do
    head -4 refcopy.fa >> temp.fa
done

#Delete spaces, newlines, tabs 
tr -d '\040\011\012\015' < temp.fa > temp.nospace.fa 

#Pick random substring of length 100 from duplicates and format it into fastq format
awk -v n=1 -v s=100  'BEGIN {"date +%N" | getline seed; srand(seed);}
                          {len=length($0);
                           for(i=1;i<=n;i++)
                              {k=rand()*(len-s)+1; printf "%s\t", substr($0,k,s)}
                               print ""}' temp.nospace.fa > temp.multi.fq
tr -d '\040\011\012\015' < temp.multi.fq > temp.multi.nospace.fq
mv temp.multi.nospace.fq temp.multi.fq
tr "[A-Z]" "\~" < temp.multi.fq > temp.multi.tildas.fq
for ((i=1;i<=3;i++));
do   
   echo "@R$i" >> multi.read.fq
   cat temp.multi.fq >> multi.read.fq
   echo "" >> multi.read.fq
   echo "+"  >> multi.read.fq
   cat temp.multi.tildas.fq >> multi.read.fq
   echo "" >>multi.read.fq
done

#Add last 2 lines of refcopy as unique lines to duplicates
tail -2 refcopy.fa >> temp.fa
#Delete spaces, newlines, tabs and attach temp as the multimapped reference file
tr -d '\040\011\012\015' < temp.fa  >> multi.ref.fa

#Add last 2 lines of refcopy as uniquely mapped read to the fastq file
tail -2 refcopy.fa > temp.multi.fq
tr -d '\040\011\012\015' < temp.multi.fq > temp.multi.nospace.fq
mv temp.multi.nospace.fq temp.multi.fq
tr "[A-Z]" "\~" < temp.multi.fq > temp.multi.tildas.fq
for ((i=1;i<=1;i++));
do
   echo "@R5" >> multi.read.fq
   cat temp.multi.fq >> multi.read.fq
   echo "" >> multi.read.fq
   echo "+"  >> multi.read.fq
   cat temp.multi.tildas.fq >> multi.read.fq
   echo "" >>multi.read.fq
done
 
#Add last 4 lines of shortReads.fq as unmapped read to fastq file
   tail -4 shortReads.fq >> multi.read.fq

 
   cp multi.ref.fa ../toy.multi.ref.fa
   cp multi.read.fq ../toy.multi.read.fq
   
rm temp*


