#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import meep as mp
import numpy as np
#import scipy as sp
#from scipy import optimize as op
from scipy import interpolate as itp
from matplotlib import pyplot as plt
from multiprocessing import Pool
# from mpl_toolkits.mplot3d import Axes3D
import meep_objects as mpo
import io
import sys
import time
#from ipywidgets import IntProgress
#from IPython.display import display

#import csv

## useful function

def convert_seconds (elapsed):
    minutes = np.floor(elapsed/60)
    secs = elapsed-minutes*60
    secs = np.round(secs*100)/100

    hours = np.int_(np.floor(minutes/60))
    minutes = np.int_(minutes-hours*60)

    return f'{hours}h-{minutes}min-{secs}s'


## function to speed up sim object creation

def sym_circular_cavity (f, df, n_back, n_groove=2, D = 0.4, DBR_period = 0.2,
                         N_rings=10, empty = False, source_pos=0, dimensions = 2, anisotropy = 0, tilt_anisotropy = 0):

    extra_space_x = 1
    domain_x = DBR_period*N_rings*2+D+extra_space_x + 1# 10*2*2*DBR_period + 1
    domain_y = domain_x
    PML = 1
    monitor_distance  = domain_x/2-extra_space_x/4

    medium_back   = mpo.anysotropic_material(n_back, anisotropy, rot_angle_3=tilt_anisotropy)
    medium_groove = mp.Medium(epsilon = n_groove**2 )

    domain = mp.Vector3(domain_x+extra_space_x+PML, domain_y+extra_space_x+PML, 0)
    device = []#mp.Block(domain, material=medium_groove)]
    # outcoupler
    # device.extend(mpo.circular_DBR_cavity(
    #     medium_back, medium_groove, D+2*DBR_period*N_rings+1,
    #     DBR_period*2, 0.5, 10,
    #     orientation = mp.Vector3(0,0,1),
    #     thickness = 0))
    # cavity
    device.extend(mpo.circular_DBR_cavity(
        medium_back, medium_groove, D,
        DBR_period, 0.5, N_rings,
        orientation = mp.Vector3(0,0,1),
        thickness = 0))
    # cavity = mpo.spiral_grating(
    #     medium_groove = medium_groove,
    #     D = D,
    #     FF = .5,
    #     DBR_period = DBR_period,
    #     N_rings = N_rings,
    #     N_arms = 0,
    #     thickness = 1,
    #     center = mp.Vector3())
    # device.extend(cavity)
    symmetries=[mp.Mirror(mp.X,phase=-1)]
    flux_or = mp.Y
    monitor_pos = mp.Vector3(0,monitor_distance,0)
    res = 100#np.int(1/(1/f/np.max([np.real(n_back),np.real(n_groove)]))*10)

    if empty:
        device = []

    pml_layers = [mp.PML(PML)]

    if df == 0 :
        source = mp.Source(mp.ContinuousSource(f,width=0.1 ),
                           component=mp.Ex,
                           center=mp.Vector3(y=source_pos) )
    else :
        source = mp.Source(mp.GaussianSource(f,df),
                           component=mp.Ex,
                           center=mp.Vector3(y=source_pos) )

    sim = mp.Simulation(cell_size=domain,
                        geometry=device,
                        sources=[source],
                        resolution=res,
                        boundary_layers=pml_layers,
                        default_material=medium_back,
                        dimensions=dimensions,
                        symmetries=symmetries,
                        eps_averaging = False)


    nfreq = 1000
    monitors = []
    fluxr = mp.FluxRegion(center=monitor_pos,
                          size=mp.Vector3(0,0,0),
                          direction=flux_or)
    monitors.append(sim.add_flux(f, df/2, nfreq, fluxr))

    if not empty :
        harminv_instance = mp.Harminv(mp.Ex, mp.Vector3(), f, df)
    else :
        harminv_instance = []

    return sim, monitors, harminv_instance


#%% function for parallel computing
def run_parallel(n_eff_h,n_eff_l,D,DBR_period, sim_end, empty=False, source_pos=0, anisotropy = 0, tilt_anisotropy = 0):

    c0 = 1
    wavelength = 0.590
    wwidth = 0.25
    f=c0/wavelength

    fmax=c0/(wavelength-wwidth/2)
    fmin=c0/(wavelength+wwidth/2)
    df=fmax-fmin

    N_periods = 30

    t0 = time.time()

    sim, monitors, H_I = sym_circular_cavity(f, df, n_back=n_eff_h,  n_groove=n_eff_l, empty=empty,
                                             D=D, DBR_period=DBR_period,
                                             N_rings = N_periods, source_pos=source_pos, dimensions = 2,
                                             anisotropy= anisotropy, tilt_anisotropy = tilt_anisotropy)

    sim.init_sim()
    fig = plt.figure(dpi=150, figsize=(10,10))
    sim.plot2D( eps_parameters={"interpolation":'none'})
    plt.show()
    fig.savefig(f'00_sourceshifted.png')

    plt.close()
    # raise Exception()
    spectra_out = []

    def print_time(sim):
        print(f'nEffH{n_eff_h}_D{int(D*1e3)}__DBRperiod{int(DBR_period*1e3)} is at {np.round(sim.round_time())} \n It has run for {convert_seconds(time.time()-t0)}\n')

    sim.run(mp.after_sources( H_I ),until=sim_end)
            # until_after_sources=mp.stop_when_fields_decayed(30, mp.Ex, mp.Vector3(), sim_end))

    for i,monitor in enumerate(monitors) :
        spectrum_f = np.array(mp.get_flux_freqs(monitor))
        spectrum = mp.get_fluxes(monitor)
        spectra_out.append(np.array(spectrum))
        plt.figure(dpi=300)
        plt.plot(1/spectrum_f, np.log10(spectra_out[i]))#
        plt.grid(True)
        plt.xlabel('wavelength')
        plt.ylabel('Transmission')
        plt.show()

    resonances_Q = []
    resonances_f = []
    resonances_amp = []
    for mode in  H_I.modes :
        if np.abs(mode.Q) > 100 :
            resonances_Q.append(np.abs(mode.Q))
            resonances_f.append(mode.freq)
            resonances_amp.append(np.abs(mode.amp))
    resonances_Q = np.array(resonances_Q)
    resonances_f = np.array(resonances_f)
    resonances_amp = np.array(resonances_amp)
    sorting = np.argsort(resonances_amp)
    resonances_Q = resonances_Q[sorting[::-1]]
    resonances_f = resonances_f[sorting[::-1]]

    t = time.time()

    print(f'\n\nThe simulation has run for {convert_seconds(t-t0)}\n')

    return [spectrum_f, spectra_out, resonances_f, resonances_Q]

#%% geometry and simulation parameters

sim_end = 300


anisotropy = 0


n_eff_l = 1
n_eff_hs = [1.1] #np.linspace(1.01,1.2,100) # [1.1]#1.0543, 1.0985, 1.1405] # 50 75 and 100 nm pmma thickness

# n_eff_h = 1
# n_eff_l = 1.1


j = 0
# for n_eff_h in n_eff_hs:
period = .280 #round(wavelength/(n_eff_l+n_eff_h),3 )
Ds = period * np.linspace(.1, 3, 100) #np.array([0.44, 2.36])#
D = .112# period * .4

# crete input vector for parallell pool. It has to be a list of tuples,
# where each element of the list represent one iteration and thus the
# element of the tuple represent the inputs.
empty = True
tuple_list = [ (n_eff_hs[0], n_eff_l,
                D, period,
                sim_end, empty,
                0,
                anisotropy,
                0 )]
empty = False


for D in Ds:
    for source_pos in [0, D/4, D/2]:
        tuple_list.append( (n_eff_hs[0], n_eff_l,
                            D, period,
                            sim_end, empty,
                            source_pos,
                            anisotropy,
                            0 ) )
    j += 1
# mp.verbosity(0)
# mp.quiet(True)
# run non parallel
output=[]
t0 = time.time()
for i in range(j):
    t1 = time.time()
    # print(tuple_list[i])
    output.append(run_parallel(*tuple_list[i]))
    print(f'It has run for {convert_seconds(time.time()-t1)}, {i+1}/{j}')
    print(f'It will take roughly {convert_seconds((time.time()-t0)/(i+1)*(j-i-1))} more')
# with Pool(5) as parfor:
#     output = parfor.starmap(run_parallel, tuple_list)
print(f'Total took {convert_seconds(time.time()-t0)}')

#%% plots
# N_resonances=0
# for var in output:
#     if len(var[2]) > N_resonances:
#         N_resonances = len(var[2])

resonance_table = []
for var in output :
    N_resonances = len(var[2])
    resonance_row = []
    for l in range(N_resonances):
        resonance_row.append([np.round(1/var[2][l]*1e3, 1), np.int(var[3][l])] )
    if N_resonances == 0 :
        resonance_row.append([ 0, 0 ])
    resonance_table.append(resonance_row)
print(resonance_table)

image = []
spectrum_empty = output[0][1][0]
for k, var in enumerate(output[1:]) :
    legend_str = []
    # resonance_table = [ [ np.round(1/var[2][l]*1e3, 1), np.int(var[3][l]) ]  for l in range(var[2].size) ]
    wavelength = 1/var[0]*1e3
    fig = plt.figure(k+1,dpi=200)
    ax = fig.add_subplot(111)
    spectrum = np.log(var[1][0])/spectrum_empty
    plt.plot(wavelength, spectrum)
    plt.xlim(wavelength.min(),wavelength.max())
    # plt.ylim(-2,2)
    ax.grid(True)
    plt.xlabel('wavelength')
    plt.ylabel('Transmission')
    plt.title(f'n_eff_h {tuple_list[k][0]}; D {int(tuple_list[k][2]*1e3)}; DBR_period {int(tuple_list[k][3]*1e3)}')
    ax2 = fig.add_subplot(336)
    # plt.title('Table of the resonances')
    collabel=[ "Wavelength", "Quality"]
    rowlabel=[ f'{i}' for i,dummy in enumerate(resonance_table[k])]#[ f'({np.int(np.rad2deg(angles[i][0]))},{np.int(np.rad2deg(angles[i][1]))})' for i in indeces]
    ax2.axis('tight')
    ax2.axis('off')
    the_table = ax2.table(cellText=resonance_table[k], colLabels=collabel, rowLabels=rowlabel,loc='center')


    plt.show()
    fig.savefig(f'nEffH{tuple_list[k][2]}_D{int(tuple_list[k][4]*1e3)}__DBRperiod{int(tuple_list[k][5]*1e3)}_sourceshifted.png')
    # plt.close(fig)
    image.append(spectrum )
image = np.array(image).transpose()
fig = mpo.plot_image(wavelength, n_eff_hs, image)
plt.xlabel('wavelength [nm]')
plt.ylabel('n_eff_h')
plt.title('Spectral response')
fig.savefig(f'n_eff_h {tuple_list[k][0]}_spacer_dependence_DBRperiod{int(tuple_list[k][3]*1e3)}_source_tilted.png')