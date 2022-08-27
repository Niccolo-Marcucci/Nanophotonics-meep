#!/bin/bash

filename="$1"
echo $filename

jobname="meep_sim"

git add $filename
git commit -m "local run: file $filename. The prefix for the simulation is the hash of this commit"

hash=$(git log -n 1 --pretty=format:"%H")
prefix="${hash:0:10}"

output_folder="${jobname}_$(date +%y%m%d-%H%M%S)_${prefix}"

prefix="$(date +%y%m%d-%H%M%S)_${prefix}"

mkdir -p data/$output_folder

mpirun -np $2 python $filename $prefix ${@:3} $2

mv *${prefix}* data/$output_folder
