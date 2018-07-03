#!/bin/bash
#Create the short read file with the correct format
rm longReads.fq*
touch longReads.fq
for ((i=1;i<=10;i++));
do
   #shuf -n 29 refcopy.fa > long1.fq
   #truncate -s 2025 long1.fq
    awk -v n=1 -v s=2000  'BEGIN {"date +%N" | getline seed; srand(seed);}
                          {len=length($0);
                           for(i=1;i<=n;i++)
                              {k=rand()*(len-s)+1; printf "%s\t", substr($0,k,s\
)}
                               print ""}' refnospace.fa > long1.fq
   #awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n"\
#);}\
#' < long1.fq > long1nospace.fq
#   mv long1nospace.fq long1.fq
   tr -d '\040\011\012\015' < long1.fq > long1shorten.fq
   mv long1shorten.fq long1.fq
   tr "[A-Z]" "\~" < long1.fq > long1_tildas.fq
   echo "@R$i" >> longReads.fq
   cat long1.fq >> longReads.fq
   echo "" >> longReads.fq
   echo "+"  >> longReads.fq
   cat long1_tildas.fq >> longReads.fq
   echo "" >> longReads.fq
done
cp longReads.fq toy.long_reads.fq
mv toy.long_reads.fq ..