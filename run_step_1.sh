#!/bin/bash

for nsims in 2 3 4 6 8 11 15 20 27 37 49 66 88 118 159 213 285 382 512

do

let x=512/nsims

for ((i=1;i<=x;i++)); 

do

matrix='minimal-'$i'-nsims'$nsims'-nclusters16-tc-lag1.00e-02'

echo $matrix

mkdir "Results/$matrix"

cp "Count_Matrices/$matrix.txt"  "Results/$matrix/"

matrix=Results/$matrix/$matrix

echo $matrix

module load matlab/2016a

matlab -r "step_1 "$matrix" 2 15; exit"

done
done
