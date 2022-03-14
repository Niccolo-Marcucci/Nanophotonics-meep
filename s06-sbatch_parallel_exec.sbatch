#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=32:00:00
# SBATCH --ntasks=64
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=1G

mpirun -np 64 python $@ 64

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
