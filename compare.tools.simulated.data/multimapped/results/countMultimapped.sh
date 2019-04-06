#!/bin/bash
for i in {10..19..1}
do
    let "numOnesinTens=numOnesinTens + $(grep "R$i" $1 | wc -l)"
done 
one=$(grep "R1" $1 | wc -l)
two=$(grep "R2" $1 | wc -l)
twenty=$(grep "R20" $1 | wc -l)
echo "$(($one - $numOnesinTens))" | tr -d '\n'
echo "," | tr -d '\n'
echo "$(($two - $twenty))" | tr -d '\n'
echo "," | tr -d '\n'
for i in {3..20..1}
do
    grep "R$i" $1 | wc -l | tr -d '\n'
    echo "," | tr -d '\n'
done
