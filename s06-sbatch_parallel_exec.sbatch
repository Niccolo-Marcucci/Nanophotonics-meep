#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=32:00:00
# SBATCH --ntasks=100
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=300M
# SBATCH --partition=global
#SBATCH --nodelist=compute-3-[13,17-20]
# the following nodes are not reserved for anyone in particular
# SBATCH --nodelist=compute-2-[5,6,9],compute-4-1

mpirun -np 100 python $@ 10

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
