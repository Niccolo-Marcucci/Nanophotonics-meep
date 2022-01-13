#!/bin/bash
#SBATCH --job-name=CW_simulaion
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=CW_simulation_%j.log
#SBATCH --mem-per-cpu=1024M

mpirun -np 24 python $1 $2

