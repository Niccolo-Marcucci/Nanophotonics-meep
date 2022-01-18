#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=00:30:00
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=200M

mpirun -np 12 python $1 $2

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
