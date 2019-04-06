#!/bin/bash
ten=$(grep "R10" $1 | wc -l)
hundred=$(grep "R100" $1 | wc -l)
echo "$(($ten - $hundred))" | tr -d '\n'
echo "," | tr -d '\n'
for i in {20..100..10}
do
    grep "R$i" $1 | wc -l | tr -d '\n'
    echo "," | tr -d '\n'
done
echo " "