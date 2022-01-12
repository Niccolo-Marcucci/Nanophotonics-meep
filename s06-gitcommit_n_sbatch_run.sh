#!/bin/bash

filename="$1"

previous_hash=$(git log -n 1 --pretty=format:"%H")
prefix="${previous_hash:0:10}"

git add $filename
git add s06-sbatch_parallel_exec.sbatch
git commit -m "run simulation $filename. The prefix for the simulation is the hash of the previous commit"

sbatch s06-sbatch_parallel_exec.sbatch $filename $prefix
