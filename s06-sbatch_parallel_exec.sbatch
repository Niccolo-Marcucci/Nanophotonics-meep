#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=32:00:00
#SBATCH --ntasks=100
# SBATCH --nodes=
# SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=300M
#SBATCH --partition=global,global*
mpirun -np 100 python $@ 100

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
