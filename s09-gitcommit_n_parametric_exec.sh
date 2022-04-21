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

njobs=$2

for (( i=0 ; i<$2 ; i++ ));
do
    python $filename ${prefix} "parallel_bash" $i $njobs 
done

wait
mv *${prefix}* data/$output_folder
