#!/bin/bash

filename="$1"

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

mkdir -p data/$prefix

jobname="${jobname}_${prefix}"
outputname="${jobname}_${prefix}_%j.log"
errorname="${jobname}_${prefix}_%j.err"

echo "$outputname, $errorname"

sbatch --job-name=$jobname --output=$outputname --error=$errorname s06-sbatch_parallel_exec.sbatch $filename $prefix
