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
import json
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


class Simulation(mp.Simulation):

    def __init__(self, sim_name='simulation_2D', dimensions=2, symmetries = []):

        self.name = sim_name

        self.extra_space_xy = .3

        self.PML_width = .6

        self._empty = True

        super().__init__(
                    cell_size = mp.Vector3(1,1,0),
                    geometry = [],
                    sources = [],
                    resolution = 1,
                    boundary_layers = [],
                    dimensions = dimensions,
                    symmetries = symmetries,
                    filename_prefix = sim_name,
                    force_complex_fields = False,
                    eps_averaging = False)

    @property
    def empty(self):
        return self._empty

    @empty.setter
    def empty(self,value):
        self._empty = value
        self.reset_meep()
        self.geometry = []
        try:
            if self._empty :
                self.geometry.extend( self._empty_geometry )
            else:
                self.geometry.extend( self._empty_geometry )
                self.geometry.extend( self._geometry )
        except AttributeError:
            raise AttributeError("cannot assign 'empty' property before initializing the geometry")

    def init_geometric_objects(self, eff_index_info={}, resolution=1, pattern_type='positive', cavity_parameters={}, outcoupler_parameters={}):
        self._geometry = []
        self._empty_geometry = []

        self.cavity_r_size = (cavity_parameters["D"]/2 + cavity_parameters["period"] * cavity_parameters["N_rings"]) * (cavity_parameters["N_rings"]>0)
        self.outcou_r_size = (outcoupler_parameters["D"]/2 + outcoupler_parameters["period"] * outcoupler_parameters["N_rings"]) * (outcoupler_parameters["N_rings"]>0)

        self.domain_x = self.domain_y = 2*(self.cavity_r_size + self.outcou_r_size + self.extra_space_xy)

        if pattern_type == 'positive':
            grating_index = np.real(eff_index_info["n_eff_l"])
            background_index = np.real(eff_index_info["n_eff_h"])
            medium_back   = mpo.anysotropic_material(background_index,
                                                     eff_index_info["anisotropy"],
                                                     rot_angle_3=eff_index_info["tilt_anisotropy"])
            medium_groove = mp.Medium(epsilon = grating_index**2 )

        elif pattern_type == 'negative':
            grating_index = np.real(eff_index_info["n_eff_h"])
            background_index = np.real(eff_index_info["n_eff_l"])
            medium_groove   = mpo.anysotropic_material(grating_index,
                                                     eff_index_info["anisotropy"],
                                                     rot_angle_3=eff_index_info["tilt_anisotropy"])
            medium_back = mp.Medium(epsilon = background_index**2 )

        else :
            raise ValueError(f'patter type "{pattern_type}" is unknown')

        self.default_material = medium_back
        if cavity_parameters["N_rings"] > 0:
            cavity = mpo.circular_DBR_cavity(
                medium_back, medium_groove,
                cavity_parameters["D"],
                cavity_parameters["period"],
                cavity_parameters["FF"],
                cavity_parameters["N_rings"],
                orientation = mp.Vector3(0,0,1),
                thickness = 0)

            self._geometry.extend(cavity)

        elif outcoupler_parameters["N_rings"] > 0:
            outcoupler = mpo.circular_DBR_cavity(
                medium_back, medium_groove,
                self.cavity_r_size*2 + outcoupler_parameters["D"],
                outcoupler_parameters["period"],
                outcoupler_parameters["FF"],
                outcoupler_parameters["N_rings"],
                orientation = mp.Vector3(0,0,1),
                thickness = 0)
            self._geometry.extend(outcoupler)

        # this  will add all geometric objects to the simulation
        self.empty = False

        # resolution is 10 points per wavelength in the highest index material time a scale factor
        self.resolution = resolution

        self.name = self.name + f'_res{self.resolution}'
        self.filename_prefix = self.name

        # round domain with an integer number of grid points
        self.grid_step = 1/self.resolution

        self.cell_size = mp.Vector3(self.domain_x + 2*self.PML_width,
                                    self.domain_y + 2*self.PML_width)
        print(self.cell_size)
        # make domain an integer number of voxels
        Nx = int(self.cell_size.x / self.grid_step)
        Nx -= np.mod(Nx,2) + 1      # make odd
        self.cell_size.x = Nx * self.grid_step
        Ny = int(self.cell_size.y / self.grid_step)
        Ny -= np.mod(Ny,2) + 1
        self.cell_size.y = Ny * self.grid_step

        print(self.cell_size)
        print()
        print(f"Number of voxels is ({Nx}x{Ny}) = {Nx*Ny/1e6} Mln")
        print(f"Minimum expected memory is {96*Nx*Ny/2**30:.2f}GB")
        print()

        self.boundary_layers = [mp.PML(self.PML_width)]
        # print( [self.cell_size.x / self.

        with open(f'{self.name}.json', 'w') as fp:
            data2save = {"eff_index_info": eff_index_info,
                         "pattern_type": pattern_type,
                         "resolution": self.resolution}

            if cavity_parameters["N_rings"] > 0:
                data2save["cavity_parameters"] = cavity_parameters

            if outcoupler_parameters["N_rings"] > 0:
                data2save["outcoupler_parameters"] = outcoupler_parameters

            json.dump(data2save, fp,  indent=4)


    def init_sources_and_monitors(self, f, df, source_pos, allow_profile=False) :
        self.sources = [ mp.Source(
            src = mp.ContinuousSource(f,fwidth=0.1) if df==0 else mp.GaussianSource(f,fwidth=df),
            center = source_pos,
            size = mp.Vector3(),
            component = mp.Ey)]

        self.harminv_instance = None
        self.field_profile = None
        self.spectrum_monitors = []

        if  allow_profile :
            self.field_profile = self.add_dft_fields([mp.Ey], f, 0, 1,
                                                     center = mp.Vector3(),
                                                     size = mp.Vector3(self.domain_x-.5*self.extra_space_xy, 0)) #, yee_grid=True))
        else:
            if self.cavity_r_size > 0 :
                DL = self.cavity_r_size + 0.02

                nfreq = 1000
                fluxr = mp.FluxRegion(
                    center = mp.Vector3(DL, 0),
                    size = mp.Vector3(0,0),
                    direction = mp.X)
                self.spectrum_monitors.append(self.add_flux(f, df, nfreq, fluxr))#, yee_grid=True))

                # if not self.empty:
                #     self.harminv_instance = mp.Harminv(mp.Ey, mp.Vector3(), f, df)

#%% function for parallel computing
def run_parallel(wavelength, n_eff_h, n_eff_l, D, DBR_period, empty=False, source_pos=0, anisotropy = 0, tilt_anisotropy = 0):
    import meep as mp

    c0 = 1
    # wavelength = 0.590
    wwidth = 0.25
    f=c0/wavelength

    sim_end=500

    fmax=c0/(wavelength-wwidth/2)
    fmin=c0/(wavelength+wwidth/2)
    df=fmax-fmin

    pattern_type = 'positive'

    t0 = time.time()

    cavity_parameters = {
        "D": D,
        "FF": .5,
        "period": DBR_period,
        "N_rings": 30}

    outcoupler_parameters = {
        "type": 'spiral',
        "D": 1,
        "FF": .5,
        "period": DBR_period * 2,
        "N_rings": 0,
        "N_arms": 0}

    eff_index_info = {
        "n_eff_h" : n_eff_h,
        "n_eff_l" : n_eff_l,
        "anisotropy" : anisotropy,
        "tilt_anisotropy" : tilt_anisotropy}


    t0 = time.time()


    date = time.strftime('%y%m%d-%H%M%S')#'211001-121139'#
    if len(sys.argv) > 1:
        sim_prefix = f"{sys.argv[1]}"
    else:
        sim_prefix = f"{date}"

    sim_name = "2D_eff_index_"
    sim_name += "cavity_" if cavity_parameters["N_rings"] > 0 else ""
    sim_name += "and_outcoupler_" if outcoupler_parameters["N_rings"] > 0 else ""
    sim_name += f"{sim_prefix}_"
    sim_name += f"D{D*1e3:.0f}_src{source_pos*1e3:.0f}"

    sim = Simulation(sim_name,symmetries=[mp.Mirror(mp.X), mp.Mirror(mp.Y,phase=-1) ])#mp.Mirror(mp.Y,phase=-1)])
    sim.extra_space_xy += wavelength/n_eff_l
    sim.eps_averaging = False
    sim.init_geometric_objects( eff_index_info = eff_index_info,
                                resolution = 100,
                                pattern_type = pattern_type,
                                cavity_parameters = cavity_parameters,
                                outcoupler_parameters = outcoupler_parameters)

    if empty:
        sim.empty = True
        sim.name += '_empty'
    else:
        sim.empty = False

    sim.init_sources_and_monitors(f, df, source_pos=mp.Vector3(x=source_pos), allow_profile=False)

    sim.init_sim()
    # fig = plt.figure(dpi=150, figsize=(10,10))
    # plot = sim.plot2D(eps_parameters={"interpolation":'none'})
    # fig.colorbar(plot.images[0])
    # # plt.show()
    # fig.savefig(f'{sim.name}-xy.jpg')
    # plt.close()
    # raise Exception()


    # mp.verbosity(0)
    sim.run(until=sim_end)
    print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to run\n')

    t = np.round(sim.round_time(), 2)

    data2save = {}
    if sim.harminv_instance != None :
        resonances_Q = []
        resonances_f = []
        for mode in  sim.harminv_instance.modes :
            if np.abs(mode.Q) > 100 :
                resonances_Q.append(np.abs(mode.Q))
                resonances_f.append(mode.freq)
        resonances_Q = np.array(resonances_Q)
        resonances_f = np.array(resonances_f)
        sorting = np.argsort(resonances_Q)
        resonances_Q = resonances_Q[sorting[::-1]]
        resonances_f = resonances_f[sorting[::-1]]

        N_resonances = len(resonances_f)
        resonance_table = []
        for l in range(N_resonances):
            resonance_table.append([np.round(1/resonances_f[l]*1e3, 1), int(resonances_Q[l])] )
        if N_resonances == 0 :
            resonance_table.append([ 0, 0 ])
        print()
        print(resonance_table)
        print()

        # with open(f'{sim.name}_output.json', 'a') as fp:
        #     data2save = {f"resonance_table_t{t}": resonance_table}
        #     json.dump(data2save, fp,  indent=4)
        data2save = {f"resonance_table_t{t}": resonance_table}

    if sim.field_profile != None:
        for j in range(sim.field_profile.nfreqs):
            data2save[f"field_profile_Ey_{j}"] = sim.get_dft_array(sim.field_profile, mp.Ey, i)
        data2save["field_profile_Eps"] = sim.get_array(mp.Dielectric,
                                                       center = sim.field_profile.regions[0].center,
                                                       size = sim.field_profile.regions[0].size)
        (x, _, _, _) = sim.get_array_metadata(center = sim.field_profile.regions[0].center,
                                              size = sim.field_profile.regions[0].size)
        data2save["field_profile_x"] = x

    spectra = []
    for monitor in sim.spectrum_monitors :
        spectrum_f = np.array(mp.get_flux_freqs(monitor))
        spectra.append(np.array(mp.get_fluxes(monitor)))

    if len(spectra) > 0 :
        data2save["wavelength"] = 1/spectrum_f*1e3
        data2save["spectra"] = spectra
    if len(data2save) > 0:
        mpo.savemat(f'{sim.name}_spectra_t{t}.mat', data2save)

    return data2save, sim.name


#%% geometry and simulation parameters
if __name__ == "__main__":              # good practise in parallel computing

    anisotropy = 0

    wavelength = .590# 0.5703#.6088#.5703#.5884#.5893#0.5947#0.5893#.5922, ]

    n_eff_l = 1
    n_eff_hs = [1.1] #np.linspace(1.01,1.2,100) # [1.1]#1.0543, 1.0985, 1.1405] # 50 75 and 100 nm pmma thickness

    # n_eff_h = 1
    # n_eff_l = 1.1

    # for n_eff_h in n_eff_hs:
    period = .280 #round(wavelength/(n_eff_l+n_eff_h),3 )
    Ds = period * np.linspace(0, 3, 500) #np.array([0, 0.45, 1, 1.5, 2.36])#0.45, 0.9, 2.36])#
    # D = .112# period * .4

    # spacers_to_test = [.08] #np.linspace(.00,.700, N)
    # output = [run_parallel("D", spacers_to_test[i], False)]

    # crete input vector for parallell pool. It has to be a list of tuples,
    # where each element of the list represent one iteration and thus the
    # element of the tuple represent the inputs.
    empty = True
    tuple_list = [ (wavelength,
                    n_eff_hs[0], n_eff_l,
                    Ds[-1], period,
                    empty,
                    0,
                    anisotropy,
                    0 )]
    empty = False

    j = 1
    for source_pos in [period/4]: # 0, period/4, period/2]:
        for D in Ds:
            tuple_list.append( (wavelength,
                                n_eff_hs[0], n_eff_l,
                                D, period,
                                empty,
                                source_pos,
                                anisotropy,
                                0 ) )
            j += 1
    mp.verbosity(0)
    # mp.quiet(True)
    # run non parallel
    output = []
    names = []
    t0 = time.time()

    try:
        from mpi4py import MPI
    except:
        non_parallel_conda = True
    else:
        non_parallel_conda = False

    if len(sys.argv) < 2 or non_parallel_conda:
        for i in range(j):
            t1 = time.time()
            # print(tuple_list[i])
            data, name = run_parallel(*tuple_list[i])
            output.append(data)
            names.append(name)
            print(f'It has run for {convert_seconds(time.time()-t1)}, {i+1}/{j}')
            print(f'It will take roughly {convert_seconds((time.time()-t0)/(i+1)*(j-i-1))} more')
            print()
            print()

    else:
        comm = MPI.COMM_WORLD
        N_jobs = int(sys.argv[-1])

        j = mp.divide_parallel_processes(N_jobs)

        N_list = len(tuple_list)
        if N_list < N_jobs :
            raise ValueError(f"Number of jobs should be lower than number of loop iterations to (f{N_list}")

        N_loops_per_job = int(N_list/N_jobs) + 1
        data_list = []
        name_list = []
        for i in range(N_loops_per_job):
            t1 = time.time()
            tuple_index = j*N_loops_per_job + i
            print(tuple_index)
            if tuple_index >= N_list :
                continue
            data, name = run_parallel(*tuple_list[tuple_index])
            data_list.append(data)
            name_list.append(name)
            print(f'It has run for {convert_seconds(time.time()-t1)}, {i+1}/{j}')
            print(f'It will take roughly {convert_seconds((time.time()-t0)/(i+1)*(j-i-1))} more')

        if mp.am_really_master():
            output.extend(data_list)
            names.extend(name_list)
            for src in range(1, N_jobs):
                output.extend( comm.recv(source=src, tag=11) )
                names.extend ( comm.recv(source=src, tag=12) )
                # comm.recv(source=src, tag=11)
                # comm.recv(source=src, tag=12)
        else:
            comm.send(data_list, dest=0, tag=11)
            comm.send(name_list, dest=0, tag=12)
            exit()
       # mp.merge_subgroup_data(output)
    # with Pool(5) as parfor:
    #     output = parfor.starmap(run_parallel, tuple_list)
    print(f'Total took {convert_seconds(time.time()-t0)}')

    #%% plots
    N_resonances=0
    # for var in output:
    #     if len(var[2]) > N_resonances:
    #         N_resonances = len(var[2])

    # resonance_table = []
    # for var in output :
    #     N_resonances = len(var['resonance_table_t300.0'])
    #     resonance_row = []
    #     for l in range(N_resonances):
    #         resonance_row.append([np.round(1/var['resonance_table_t300.0'][l][0]*1e3, 1), np.int(var['resonance_table_t300.0'][l][1])] )
    #     if N_resonances == 0 :
    #         resonance_row.append([ 0, 0 ])
    #     resonance_table.append(resonance_row)
    # print(resonance_table)
    image = []
    spectrum_empty = output[0]["spectra"][0]
    l=0
    for k, var in enumerate(output[1:]) :
        k=k+1
        spectrum = ( var['spectra'][0]/spectrum_empty)
        wavelength = var["wavelength"]
        # l=np.int((k-1)/3)
        # # resonance_table = [ [ np.round(1/var[2][l]*1e3, 1), np.int(var[3][l]) ]  for l in range(var[2].size) ]
        # if np.mod(k-1,3) == 0 :
        #     fig = plt.figure(l,dpi=150,figsize=(10,5))
        #     ax = fig.add_subplot(111)
        #     legend_str = []
        # ax.plot(wavelength, spectrum)
        # plt.xlim(550,650)
        # # plt.ylim(-2,2)
        # ax.grid(True)
        # plt.xlabel('Wavelength [nm]')
        # plt.ylabel('Transmission')
        # legend_str.append(f"sourcePos {tuple_list[k][6]*1e3:.0f}")
        # plt.title(f'n_eff_h={tuple_list[k][1]:.2f};   DBR_period={tuple_list[k][4]*1e3:.0f};   D={tuple_list[k][3]/tuple_list[k][2]:.2f}*DBR_period')
        # # ax2 = fig.add_subplot(336)
        # # # plt.title('Table of the resonances')
        # # collabel=[ "Wavelength", "Quality"]
        # # rowlabel=[ f'{i}' for i,_ in enumerate(var["resonance_table_t300.0"])]#[ f'({np.int(np.rad2deg(angles[i][0]))},{np.int(np.rad2deg(angles[i][1]))})' for i in indeces]
        # # ax2.axis('tight')
        # # ax2.axis('off')
        # # the_table = ax2.table(cellText=var["resonance_table_t300.0"], colLabels=collabel, rowLabels=rowlabel,loc='center')
        # plt.legend(legend_str)

        # # plt.show()
        # fig.savefig(f'{names[k]}_spectrum.png')
        # plt.close(fig)
        image.append(np.log10(spectrum))

        # fig = plt.figure(10*k,dpi=150,figsize=(10,5))
        # ax = fig.add_subplot(111)
        # ax.plot(var["field_profile_x"], var["field_profile_Eps"])
        # field = var["field_profile_Ey_0"]
        # field = np.abs(field/np.max(field))**2
        # ax.plot(var["field_profile_x"], field)
        # ax.grid(True)
        # plt.xlabel('x [um]')
        # plt.legend(["dielectric_constant", "field_profile"])
        # plt.title(f'n_eff_h={tuple_list[k][1]:.2f};   DBR_period={tuple_list[k][4]*1e3:.0f}; D={tuple_list[k][3]/tuple_list[k][4]:.2f}*DBR_period')
        # fig.savefig(f'{names[k]}_spectrumfield_profile.png')


    image = np.array(image).transpose()
    fig = mpo.plot_image(wavelength*1e-3, Ds, image)
    fig.set_figheight(6)
    fig.set_figwidth(12)
    fig.set_dpi(150)
    plt.xlabel('wavelength [um]')
    plt.ylabel('Spacer [um]')
    plt.title('Spectral response')
    fig.savefig(f'{names[0]}_spacer_dependence_DBRperiod{period*1e3:.0f}_sourcePos{source_pos}.png')
