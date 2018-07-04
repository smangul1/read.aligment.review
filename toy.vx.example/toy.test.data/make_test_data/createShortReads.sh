#!/bin/bash
#Create short read files with the correct format

#Delete old short read files
rm shortReads.fq*
touch shortReads.fq
for ((i=1;i<=10;i++));
do
#Choose random 100 consecutive bp from refnospace (reference file without fasta header)
  awk -v n=1 -v s=100  'BEGIN {"date +%N" | getline seed; srand(seed);} 
                          {len=length($0); 
                           for(i=1;i<=n;i++) 
                              {k=rand()*(len-s)+1; printf "%s\t", substr($0,k,s)}
                               print ""}' refnospace.fa > short1.fq   
# Remove spaces, tabs, newlines
   tr -d '\040\011\012\015' < short1.fq > short1shorten.fq
   mv short1shorten.fq short1.fq
# Create a read in a fasta format
   echo ">R$i" >> shortReads.fa
   cat short1.fq >> shortReads.fa
   echo "" >> shortReads.fa
# Create a read in a fastq format
   tr "[A-Z]" "\~" < short1.fq > short1_tildas.fq
   echo "@R$i" >> shortReads.fq
   cat short1.fq >> shortReads.fq
   echo "" >> shortReads.fq
   echo "+"  >> shortReads.fq
   cat short1_tildas.fq >> shortReads.fq
   echo "" >> shortReads.fq
done
cp shortReads.fq toy.short_reads.fq
cp shortReads.fa toy.short_reads.fa
mv toy.short_reads.fq ..
mv toy.short_reads.fa ..
cp reference.fa toy.reference.fa
mv toy.reference.fa ..
