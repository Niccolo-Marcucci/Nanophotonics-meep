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

outputname="${jobname}${empty}_${prefix}_%j.log"
errorname="${jobname}${empty}_${prefix}_%j.err"

jobname="${jobname}${empty}_${prefix}"
prefix="$(date +%y%m%d-%H%M%S)_${prefix}"

mkdir -p data/$output_folder

echo "$outputname, $errorname"

njobs=$2

for (( i=0 ; i<$2 ; i++ ));
do
    sbatch --job-name=$jobname --output=$outputname --error=$errorname s09-sbatch_single_exec.sbatch $filename $prefix "parallel_bash" $i $njobs
done