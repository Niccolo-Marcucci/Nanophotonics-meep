#!/bin/bash
#SBATCH --job-name=CW_simulaion
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=05:00:00
#SBATCH --nodes=1
# #SBATCH --ntasks=40
# #SBATCH --cpus-per-task=1
#SBATCH --output=CW_simulation_%j.log
#SBATCH --error=CW_simulation_%j.err
#SBATCH --mem=128G

mpirun -np 20 python $1 $2

