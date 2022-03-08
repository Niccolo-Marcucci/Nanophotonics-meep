#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=12:00:00
#SBATCH --ntasks=32
# SBATCH --nodes=1
# SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=500M

mpirun -np 32 python $1 32 ${@:2}

mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
