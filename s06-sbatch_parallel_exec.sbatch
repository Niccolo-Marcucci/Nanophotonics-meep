#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=07-00:00:00
#SBATCH --ntasks=64
# SBATCH --nodes=2
# SBATCH --ntasks-per-node=
#SBATCH --mem-per-cpu=200M
# SBATCH --partition=global
# SBATCH --nodelist=compute-3-[13,17-20]
# the following nodes are not reserved for anyone in particular
# SBATCH --nodelist=compute-2-[5,6,9],compute-4-1

mpirun -np 64 python $@ 64

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2

rm $1
