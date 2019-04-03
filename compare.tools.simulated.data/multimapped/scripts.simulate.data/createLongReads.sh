#!/bin/bash
#Create the short read file with the correct format
rm longReads.fq*
rm longReads.fa*
touch longReads.fq
for ((i=1;i<=10;i++));
do
    awk -v n=1 -v s=2000  'BEGIN {"date +%N" | getline seed; srand(seed);}
                          {len=length($0);
                           for(i=1;i<=n;i++)
                              {k=rand()*(len-s)+1; printf "%s\t", substr($0,k,s\
)}
                               print ""}' refnospace.fa > long1.fq
   tr -d '\040\011\012\015' < long1.fq > long1shorten.fq
   mv long1shorten.fq long1.fq
   echo ">R$i" >> longReads.fa
   cat long1.fq >> longReads.fa
   echo "" >> longReads.fa
   tr "[A-Z]" "\~" < long1.fq > long1_tildas.fq
   echo "@R$i" >> longReads.fq
   cat long1.fq >> longReads.fq
   echo "" >> longReads.fq
   echo "+"  >> longReads.fq
   cat long1_tildas.fq >> longReads.fq
   echo "" >> longReads.fq
done

rm long1*