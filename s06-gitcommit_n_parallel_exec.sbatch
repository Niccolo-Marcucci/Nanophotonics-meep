#!/bin/bash
#SBATCH --job-name=meep_run_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$2
#SBATCH --output=meep_run_test_%j.log
#SBATCH --mem-per-cpu=512M

filename="$1"

previous_hash=$(git log -n 1 --pretty=format:"%H")
prefix="${previous_hash:0:10}"

git add $filename
git commit -m "run simulation $filename. The prefix for the simulation is the hash of the previous commit"

mpirun -np $2 python $filename $prefix
