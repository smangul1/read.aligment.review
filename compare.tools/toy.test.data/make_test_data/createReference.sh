#!/bin/bash
rm reference.fa
rm refcopy.fa
#create random file - 256 byte values to 64 characters
dd if=/dev/urandom bs=8000 count=1 | base64 > random.fa
tr "[A-P]" "A"  < random.fa > refA.fa
tr "[Q-Za-f]" "C"  < refA.fa > refC.fa
tr "[g-v]" "T"  < refC.fa > refT.fa
tr "[w-z0-9+=]" "G"  < refT.fa > refG.fa
tr -d '/' < refG.fa > refcopy.fa
truncate -s 10375 refcopy.fa #truncate to 10kb
#Delete newlines so sequence is in one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < refcopy.fa > refnospace.fa
touch reference.fa
echo ">Reference" >> reference.fa
cat refnospace.fa >> reference.fa
rm refA.fa
rm refC.fa
rm refT.fa
rm refG.fa
rm random.fa