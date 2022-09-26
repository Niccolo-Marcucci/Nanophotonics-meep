#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=niccolo.marcucci@polito.it
#SBATCH --partition=global
#SBATCH --time=00-07:00:00
#SBATCH --ntasks=48
# SBATCH --nodes=2
# SBATCH --ntasks-per-node=
#SBATCH --mem-per-cpu=500M
# SBATCH --partition=global
# SBATCH --nodelist=compute-3-[13,17-20]
# the following nodes are not reserved for anyone in particular
# SBATCH --nodelist=compute-2-[5,6,9],compute-4-1

mpirun -np 48 python $@ 48

wait

python mv_files_to_folder_parallelly.py $2 48

# mv -v $(ls *${2}* | grep -v 'log\|err') data/*$2
