#!/bin/bash

filename="$1"
echo $filename

previous_hash=$(git log -n 1 --pretty=format:"%H")
prefix="${previous_hash:0:10}"

git add $filename
git commit -m "run simulation $filename. \nthe prefix for the simulation is the hash of the previous commit"

mpirun -np $2 python $filename $prefix
