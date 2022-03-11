#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=24:00:00
#SBATCH --ntasks=84
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=1500M

mpirun -np 84 python $@ 84

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
