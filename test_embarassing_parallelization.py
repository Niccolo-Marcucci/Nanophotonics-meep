#!/usr/bin/env python3

import meep as mp
import numpy as np
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
import sys

#%% function for parallel computing
def run_parallel(x, y):

    cell = mp.Vector3(16,16,0)
    geometry = [mp.Block(mp.Vector3(12,1,mp.inf),
                         center=mp.Vector3(-2.5,-3.5),
                         material=mp.Medium(epsilon=12)),
                mp.Block(mp.Vector3(1,12,mp.inf),
                         center=mp.Vector3(3.5,2),
                         material=mp.Medium(epsilon=12))]
    pml_layers = [mp.PML(1.0)]


    sources = [mp.Source(mp.ContinuousSource(wavelength=2*(11**0.5), width=20),
                         component=mp.Ez,
                         center=mp.Vector3(-7,-3.5),
                         size=mp.Vector3(0,1))]

    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        geometry=geometry,
                        sources=sources,
                        resolution=20)

    sim.run(until=100)

def proxy_fun( tuple_ ):
    run_parallel( *tuple_ )

#%% geometry and simulation parameters
if __name__ == "__main__":              # good practise in parallel computing

    x = np.linspace(0,1,6)
    y = 0
    tuple_list = [ (xx, y) for xx in x]

    N_jobs = mp.count_processors()

    if len(sys.argv) > 1 :
        if sys.argv[1] == "serial" :
            for i in range(len(tuple_list)):
                run_parallel(*tuple_list[i])

        elif sys.argv[1] == "bash" :
            j = int(sys.argv[2])

            N_list = len(tuple_list)
            N_jobs = 6
            N_loops_per_job = int(np.ceil(N_list/N_jobs))

            for i in range(N_loops_per_job):
                tuple_index = N_loops_per_job*i + j
                print(tuple_index)
                if tuple_index >= N_list :
                    continue
                else:
                    run_parallel(*tuple_list[tuple_index])

    elif N_jobs == 1 :
        with ProcessPoolExecutor(6) as parfor:
          parfor.map(proxy_fun, tuple_list)

    else:
        # parallel execution using mpirun
        j = mp.divide_parallel_processes(N_jobs)

        N_list = len(tuple_list)

        N_loops_per_job = int(np.ceil(N_list/N_jobs))

        for i in range(N_loops_per_job):
            tuple_index = N_loops_per_job*i + j
            print(tuple_index)
            if tuple_index >= N_list :
                continue
            else:
                run_parallel(*tuple_list[tuple_index])

