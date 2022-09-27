#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from mpi4py import MPI

prefix = sys.argv[1] # "lkn1wlm192131"

files = os.listdir("./data")
for file in files :
    if file.find( prefix ) >= 0 and os.path.isdir('data/' + file):
        folder = 'data/' + file

files = os.listdir("./")
def run_parallel(file):
    if file.find(prefix) >= 0:
        if file[-4:] != ".log" and file[-4:] != ".err" :
            os.system(f"mv -v ./{file} --target-directory={folder}")


# mp.reset_meep()
comm = MPI.COMM_WORLD
N_jobs = int(sys.argv[-1])

j = comm.Get_rank()

N_list = len(files)
if N_list < N_jobs :
    raise ValueError(f"Number of jobs should be lower than number of loop iterations ({N_list})")

reminder = np.mod(N_list,N_jobs)
N_loops_per_job = int(N_list/N_jobs)
if j < reminder:
    N_loops_per_job += 1

data_list = []
name_list = []
for i in range(N_loops_per_job):
    if j < reminder:
        index = j*N_loops_per_job + i
    else:
        index = reminder*(N_loops_per_job+1) + (j-reminder)*N_loops_per_job + i

    if index >= N_list :
        continue

    run_parallel(files[index])
    if j == 0:
        print(i)

