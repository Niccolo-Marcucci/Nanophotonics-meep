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
                    cell_size = mp.Vector3(1,1,1),
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


    def init_sources_and_monitors(self, f, df, source_pos) :
        self.sources = [ mp.Source(
            src = mp.ContinuousSource(f,fwidth=0.1) if df==0 else mp.GaussianSource(f,fwidth=df),
            center = source_pos,
            size = mp.Vector3(),
            component = mp.Ey)]

        self.harminv_instance = None
        self.spectrum_monitors = []

        if self.cavity_r_size > 0 :
            DL = self.cavity_r_size + 0.02

            nfreq = 1000
            fluxr = mp.FluxRegion(
                center = mp.Vector3(DL, 0),
                size = mp.Vector3(0,0),
                direction = mp.X)
            self.spectrum_monitors.append(self.add_flux(f, df, nfreq, fluxr))#, yee_grid=True))

            if not self.empty:
                self.harminv_instance = mp.Harminv(mp.Ey, mp.Vector3(), f, df)

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
def run_parallel(n_eff_h,n_eff_l,D,DBR_period, empty=False, source_pos=0, anisotropy = 0, tilt_anisotropy = 0):
    import meep as mp

    c0 = 1
    wavelength = 0.590
    wwidth = 0.25
    f=c0/wavelength

    sim_end=300

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
        "tilt_anisotropy" : tilt_anisotropy,}


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
    sim_name += f"D_{D*1e3:.0f}_source_{source_pos*1e3:.0f}"



    sim = Simulation(sim_name,symmetries=[])#mp.Mirror(mp.Y,phase=-1)])
    sim.extra_space_xy += wavelength/n_eff_l
    sim.eps_averaging = False
    sim.init_geometric_objects( eff_index_info = eff_index_info,
                                resolution = 50,
                                pattern_type = pattern_type,
                                cavity_parameters = cavity_parameters,
                                outcoupler_parameters = outcoupler_parameters)


    if empty:
        sim.empty = True
        sim.name += '_empty'
    else:
        sim.empty = False


    sim.init_sources_and_monitors(f, df, source_pos=mp.Vector3(x=source_pos))

    sim.init_sim()
    fig = plt.figure(dpi=150, figsize=(10,10))
    plot =  sim.plot2D( eps_parameters={"interpolation":'none'})
    fig.colorbar(plot.images[0])
    # plt.show()
    fig.savefig(f'{sim.name}-xy.jpg')
    plt.close()
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

        with open(f'{sim.name}_output.json', 'a') as fp:
            data2save = {f"resonance_table_t{t}": resonance_table}
            json.dump(data2save, fp,  indent=4)

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


    n_eff_l = 1
    n_eff_hs = [1.1] #np.linspace(1.01,1.2,100) # [1.1]#1.0543, 1.0985, 1.1405] # 50 75 and 100 nm pmma thickness

    # n_eff_h = 1
    # n_eff_l = 1.1

    # for n_eff_h in n_eff_hs:
    period = .280 #round(wavelength/(n_eff_l+n_eff_h),3 )
    Ds = period *np.array([0.45, 0.9, 2.36])#np.linspace(.1, 3, 4) #
    # D = .112# period * .4

    # spacers_to_test = [.08] #np.linspace(.00,.700, N)
    # output = [run_parallel("D", spacers_to_test[i], False)]

    # crete input vector for parallell pool. It has to be a list of tuples,
    # where each element of the list represent one iteration and thus the
    # element of the tuple represent the inputs.
    empty = True
    tuple_list = [ (n_eff_hs[0], n_eff_l,
                    Ds[-1], period,
                    empty,
                    0,
                    anisotropy,
                    0 )]
    empty = False

    j = 1
    for D in Ds:
        for source_pos in [0, D/4, D/2]:
            tuple_list.append( (n_eff_hs[0], n_eff_l,
                                D, period,
                                empty,
                                source_pos,
                                anisotropy,
                                0 ) )
            j += 1
    # mp.verbosity(0)
    # mp.quiet(True)
    # run non parallel
    output = []
    names = []
    t0 = time.time()

    # try:
    #     from mpi4py import MPI
    # except:
    for i in range(j):
        t1 = time.time()
        # print(tuple_list[i])
        data, name = run_parallel(*tuple_list[i])
        output.append(data)
        names.append(name)
        print(f'It has run for {convert_seconds(time.time()-t1)}, {i+1}/{j}')
        print(f'It will take roughly {convert_seconds((time.time()-t0)/(i+1)*(j-i-1))} more')
    # else:
    #     comm = MPI.COMM_WORLD
    #     i = mp.divide_parallel_processes(int(sys.argv[1]))
    #     data, name = run_parallel(*tuple_list[i])
    #     if mp.am_really_master():
    #         output.append(data)
    #         names.append(name)
    #         for src in range(1,int(sys.argv[1])):
    #             output.append( comm.recv(source=src,tag=11) )
    #             names.append ( comm.recv(source=src,tag=12) )
    #     else:
    #         comm.send(data,tag=11)
    #         comm.send(data,tag=12)
    #         exit()
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
    for k, var in enumerate(output[1:]) :
        k=k+1
        legend_str = []
        # resonance_table = [ [ np.round(1/var[2][l]*1e3, 1), np.int(var[3][l]) ]  for l in range(var[2].size) ]
        wavelength = var["wavelength"]
        fig = plt.figure(k,dpi=200)
        ax = fig.add_subplot(111)
        spectrum = ( var['spectra'][0]/spectrum_empty)
        plt.plot(wavelength, spectrum)
        plt.xlim(wavelength.min(),wavelength.max())
        # plt.ylim(-2,2)
        ax.grid(True)
        plt.xlabel('wavelength')
        plt.ylabel('Transmission')
        plt.title(f'n_eff_h {tuple_list[k][0]:.2f}; D {tuple_list[k][2]*1e3:.0f}; DBR_period {tuple_list[k][3]*1e3:.0f}, sourcePos {tuple_list[k][5]*1e3:.0f}')
        ax2 = fig.add_subplot(336)
        # plt.title('Table of the resonances')
        collabel=[ "Wavelength", "Quality"]
        rowlabel=[ f'{i}' for i,_ in enumerate(var["resonance_table_t300.0"])]#[ f'({np.int(np.rad2deg(angles[i][0]))},{np.int(np.rad2deg(angles[i][1]))})' for i in indeces]
        ax2.axis('tight')
        ax2.axis('off')
        the_table = ax2.table(cellText=var["resonance_table_t300.0"], colLabels=collabel, rowLabels=rowlabel,loc='center')


        plt.show()
        fig.savefig(f'names[k]_spectrum.png')
        # plt.close(fig)
        image.append(spectrum )
    # image = np.array(image).transpose()
    # fig = mpo.plot_image(wavelength, n_eff_hs, image)
    # plt.xlabel('wavelength [nm]')
    # plt.ylabel('n_eff_h')
    # plt.title('Spectral response')
    # fig.savefig(f'n_eff_h {tuple_list[k][0]}_spacer_dependence_DBRperiod{int(tuple_list[k][3]*1e3)}_source_tilted.png')