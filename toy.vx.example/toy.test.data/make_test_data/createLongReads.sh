#!/bin/bash
#Create short read files in the correct format

#delete old longReads 
rm longReads.fq*
rm longReads.fa*
#create a blank longReads.fq 
touch longReads.fq

for ((i=1;i<=10;i++));
do
#Choose random 2000 consecutive bp from refnospace (reference file without fasta header)
    awk -v n=1 -v s=2000  'BEGIN {"date +%N" | getline seed; srand(seed);}
                          {len=length($0);
                           for(i=1;i<=n;i++)
                              {k=rand()*(len-s)+1; printf "%s\t", substr($0,k,s\
)}
                               print ""}' refnospace.fa > long1.fq
# Remove spaces, tabs, newlines
   tr -d '\040\011\012\015' < long1.fq > long1shorten.fq
   mv long1shorten.fq long1.fq
# Create a read in a fasta format
   echo ">R$i" >> longReads.fa
   cat long1.fq >> longReads.fa
   echo "" >> longReads.fa
# Create a read in a fastq format
   tr "[A-Z]" "\~" < long1.fq > long1_tildas.fq
   echo "@R$i" >> longReads.fq
   cat long1.fq >> longReads.fq
   echo "" >> longReads.fq
   echo "+"  >> longReads.fq
   cat long1_tildas.fq >> longReads.fq
   echo "" >> longReads.fq
done
cp longReads.fq toy.long_reads.fq
cp longReads.fa toy.long_reads.fa
mv toy.long_reads.fa ..
mv toy.long_reads.fq ..
