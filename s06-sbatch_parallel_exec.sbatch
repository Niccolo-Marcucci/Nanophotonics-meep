#!/bin/bash
#SBATCH --job-name=meep_run_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --output=meep_run_test_%j.log
#SBATCH --mem-per-cpu=2048M

mpirun -np 12 python $1 $2

