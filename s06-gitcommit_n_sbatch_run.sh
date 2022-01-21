#!/bin/bash

filename="$1"

if [ $# -gt 2 ]
then
    empty="_$3"
else
    empty=""
fi
if [ $# -gt 1 ]
then
    jobname="$2"
else
    jobname="meep_sim"
fi


git add $filename
git add s06-sbatch_parallel_exec.sbatch
git commit -m "cluster run: job $jobname, file $filename. The prefix for the simulation is the hash of this commit"

hash=$(git log -n 1 --pretty=format:"%H")
prefix="${hash:0:10}"

output_folder="${jobname}_$(date +%y%m%d-%H%M%S)_${prefix}"

outputname="${jobname}_${prefix}${empty}_%j.log"
errorname="${jobname}_${prefix}${empty}_%j.err"

jobname="${jobname}_${prefix}${empty}"
prefix="$(date +%y%m%d-%H%M%S)_${prefix}"

mkdir -p data/$output_folder

echo "$outputname, $errorname"

sbatch --job-name=$jobname --output=$outputname --error=$errorname s06-sbatch_parallel_exec.sbatch $filename $prefix ${@:3}
