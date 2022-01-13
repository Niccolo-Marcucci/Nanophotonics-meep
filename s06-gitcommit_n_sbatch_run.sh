#!/bin/bash

filename="$1"

if [ $# -gt 1 ]
then
    jobname="$2"
else
    jobname="meep_sim"
fi

previous_hash=$(git log -n 1 --pretty=format:"%H")
prefix="${previous_hash:0:10}"

git add $filename
git add s06-sbatch_parallel_exec.sbatch
git commit -m "cluster run: job $jobname, file $filename, prefix $prefix. The prefix for the simulation is the hash of the previous commit"

sbatch --job-name=$jobname --output=$jobname_$prefix_%j.log --error=$jobname_$prefix_%j.err s06-sbatch_parallel_exec.sbatch $filename $prefix
