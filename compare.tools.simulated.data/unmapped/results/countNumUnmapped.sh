#!/bin/bash
rm *.junc
rm *.temp
rm *.vcf
unmapped="$(grep -l 'unmapped' short* | wc -l)" #total unmapped
total="$(ls short* | wc -l)" #total files
noUnmapped="$(grep -L 'unmapped' short* | wc -l)"
echo "unmapped = $unmapped"
grep -l 'unmapped' short* 
echo "no unmapped = $noUnmapped"
grep -L 'unmapped' short* 
echo "total = $total"
numShort="$(grep -l 'simple_short' short* | wc -l)" #number of tools that can read short reads
echo "short = $numShort"
grep -l 'simple_short' short*
numLong="$(grep -l 'simple_long' short* | wc -l)" 
echo "long = $numLong"
grep -l 'simple_long' short*