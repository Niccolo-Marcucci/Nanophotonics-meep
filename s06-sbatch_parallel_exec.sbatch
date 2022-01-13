#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=04:00:00
#SBATCH --ntasks=32
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --mem=35G

mpirun -np 32 python $1 $2

