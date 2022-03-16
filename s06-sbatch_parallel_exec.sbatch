#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=32:00:00
# SBATCH --ntasks=128
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=300M

mpirun -np 96 python $@ 96

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
