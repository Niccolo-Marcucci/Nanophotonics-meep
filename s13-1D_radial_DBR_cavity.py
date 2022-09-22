#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import meep as mp
import numpy as np
#import scipy as sp
#from scipy import optimize as op
from scipy import interpolate as itp, io
from matplotlib import pyplot as plt
from multiprocessing import Pool
# from mpl_toolkits.mplot3d import Axes3D
import meep_objects as mpo
import json
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


    def __init__(self, sim_name='simulation_2D', dimensions=mp.CYLINDRICAL, symmetries = []):

        self.name = sim_name

        self.extra_space_xy = .5

        self.PML_width = .6

        self._empty = True

        self.epsilon_proxy_function = lambda pos: self.circular_undeformed_cavity(pos) #imported_structure(pos) #

        super().__init__(
                    cell_size = mp.Vector3(1,1,0),
                    geometry = [],
                    sources = [],
                    resolution = 1,
                    boundary_layers = [],
                    dimensions = dimensions,
                    m = 0,
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
                self.epsilon_func = None
            else:
                self.geometry.extend( self._empty_geometry )
                self.epsilon_func = self.epsilon_proxy_function
        except AttributeError:
            raise AttributeError("cannot assign 'empty' property before initializing the geometry")

    def init_geometric_objects(self, eff_index_info={}, resolution=1, pattern_type='positive', cavity_parameters={}):
        self._empty_geometry = []
        self.cavity_parameters = cavity_parameters
        self.eff_index_info = eff_index_info
        self.cavity_r_size = (cavity_parameters["D"]/2 + cavity_parameters["period"] * cavity_parameters["N_rings"]) * (cavity_parameters["N_rings"]>0)

        self.domain_x = self.cavity_r_size + self.extra_space_xy

        if pattern_type == 'positive':
            self.grating_index = np.real(eff_index_info["n_eff_l"])
            self.background_index = np.real(eff_index_info["n_eff_h"])
            self.medium_back   = mpo.anisotropic_material(self.background_index,
                                                     eff_index_info["anisotropy"],
                                                     rot_angle_3=eff_index_info["tilt_anisotropy"])
            self.medium_groove = mpo.Medium(epsilon = self.grating_index**2 )

        elif pattern_type == 'negative':
            self.grating_index = np.real(eff_index_info["n_eff_h"])
            self.background_index = np.real(eff_index_info["n_eff_l"])
            self.medium_groove   = mpo.anisotropic_material(self.grating_index,
                                                     eff_index_info["anisotropy"],
                                                     rot_angle_3=eff_index_info["tilt_anisotropy"])
            self.medium_back = mpo.Medium(epsilon = self.background_index**2 )

        else :
            raise ValueError(f'patter type "{pattern_type}" is unknown')

        self.default_material = self.medium_back

        self.epsilon_func = self.epsilon_proxy_function # this assignement is probably unnecessary due to the empty variable

        # this  will add all geometric objects to the simulation
        self.empty = False

        # resolution is 10 points per wavelength in the highest index material time a scale factor
        self.resolution = resolution

        self.name = self.name + f'_res{self.resolution}'
        self.filename_prefix = self.name

        # round domain with an integer number of grid points
        self.grid_step = 1/self.resolution

        self.cell_size = mp.Vector3(self.domain_x + self.PML_width)
        print(self.cell_size)
        # make domain an integer number of voxels
        Nx = int(self.cell_size.x / self.grid_step)
        Nx -= np.mod(Nx,2) # make even; + 1      # make odd
        self.cell_size.x = Nx * self.grid_step

        print(self.cell_size)
        print()
        print(f"Number of voxels is ({Nx}) = {Nx/1e6} Mln")
        print(f"Minimum expected memory is {96*Nx/2**30:.3f}GB")
        print()

        self.boundary_layers = [mp.PML(self.PML_width)]
        # print( [self.cell_size.x / self.

        with open(f'{self.name}.json', 'w') as fp:
            data2save = {"eff_index_info": {key: eff_index_info[key] for key in eff_index_info.keys() if type(eff_index_info[key]) != type(lambda x:x) },
                         "pattern_type": pattern_type,
                         "resolution": self.resolution}

            if cavity_parameters["N_rings"] > 0:
                data2save["cavity_parameters"] = cavity_parameters

            json.dump(data2save, fp,  indent=4)

    def circular_undeformed_cavity(self, pos):
        r = pos.x
        D = self.cavity_parameters["D"]
        FF = self.cavity_parameters["FF"]
        period = self.cavity_parameters["period"]
        N = self.cavity_parameters["N_rings"]

        is_groove = False
        for i in range(N):
            groove_start = D/2 + i*period
            groove_end   = D/2 + FF*period + i*period
            if r > groove_start and r <= groove_end:
                is_groove = True
                break

        if is_groove:
            local_index = self.grating_index # + mod_tranches
        else:
            local_index = self.background_index # + mod_ridges

        if r < D/2 or r > D/2 + N*period - (1-FF)*period:
            local_index = self.eff_index_info["spacer_index"]

            # Z = self.weird_cone( pos)
            # local_index = np.polyval(p_neff_590, Z+60.8)

        return local_index**2

    def init_sources_and_monitors(self, f, df, source_pos, source_tilt, allow_profile=False) :
        self.sources = [ mp.Source(
                            src = mp.ContinuousSource(f,fwidth=0.1) if df==0 else mp.GaussianSource(f,fwidth=df),#,,is_integrated=True
                            center = source_pos,
                            size = mp.Vector3(y = 0), #self.cell_size.y),#
                            component = mp.Ez,
                            amplitude = np.cos(source_tilt))] #

        self.harminv_instance = None
        self.field_profile = None
        self.field_FT = None
        self.spectrum_monitors = []
        self.time_monitors = []
        self.Ex = []
        self.Ey = []
        self.Ez = []

        if  allow_profile :
            self.field_profile = self.add_dft_fields([mp.Ez], 1/np.array([.5772, .5842, .5854]),#f, 0, 1,
                                                     center = mp.Vector3(),
                                                     size = mp.Vector3(self.domain_x-.5*self.extra_space_xy)) #, yee_grid=True))
        else:
            if self.cavity_r_size > 0 :
                DL = self.cavity_r_size + self.extra_space_xy*.5
                nfreq = 1000 if df != 0 else 1
                DL_x = DL
                direction = mp.R
                fluxr = mp.FluxRegion(
                    center = mp.Vector3(DL_x),
                    size = mp.Vector3(0, 0),
                    direction = direction)
                self.spectrum_monitors.append(self.add_flux(f, df, nfreq, fluxr))#, yee_grid=True))
                # self.time_monitors.append(mp.Volume(center = mp.Vector3(DL_x, DL_y), size = mp.Vector3(0, 0)))
                self.Ex.append([])
                self.Ey.append([])
                self.Ez.append([])
                # self.field_FT = self.add_dft_fields([mp.Ez], f, df, nfreq,
                #                                     center = mp.Vector3(self.cavity_parameters["D"]/2),
                #                                     size = mp.Vector3(self.cavity_parameters["D"]))#self.cavity_parameters["D"]/2,self.cavity_parameters["D"]/2 ))
                self.time_monitors.append(mp.Volume(center = mp.Vector3(),
                                                    size = mp.Vector3(0)))
                self.Ex.append([])
                self.Ey.append([])
                self.Ez.append([])

                if not self.empty:
                    self.harminv_instance = None # mp.Harminv(mp.Ez, mp.Vector3(), f, df)

def save_fields(sim):
    i=-1
    for i, monitor in enumerate(sim.time_monitors):
        sim.Ex[i].append( sim.get_array(mp.Ex, center = monitor.center, size = monitor.size) )
        sim.Ey[i].append( sim.get_array(mp.Ey, center = monitor.center, size = monitor.size) )
        sim.Ez[i].append( sim.get_array(mp.Ez, center = monitor.center, size = monitor.size) )

#%% function for parallel computing
def run_parallel(wavelength, n_eff_h, n_eff_l, n_eff_spacer, D, DBR_period, empty=False, source_pos=0, source_tilt=0, n_eff_mod_l = 0, n_eff_mod_h = 0, n_eff_wv=None, Z_f=None):
    # import meep as mp

    c0 = 1
    # wavelength = 0.590
    wwidth = 0.15
    f=c0/wavelength

    sim_end=400

    fmax=c0/(wavelength-wwidth/2)
    fmin=c0/(wavelength+wwidth/2)
    df=fmax-fmin

    pattern_type = 'positive'

    t0 = time.time()

    cavity_parameters = {
        "D": D,
        "FF": .48,
        "period": DBR_period,
        "N_rings": 30,
        "tilt": 0} #source_tilt}

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
        "anisotropy" : 0,
        "tilt_anisotropy" : 0,
        "modulation_amplitude_ridges": n_eff_mod_h,
        "modulation_amplitude_tranches": n_eff_mod_l,
        "spacer_index": n_eff_spacer,
        "n_eff_wv": n_eff_wv,
        "Z_f": Z_f}


    t0 = time.time()


    date = time.strftime('%y%m%d-%H%M%S')#'211001-121139'#
    if len(sys.argv) > 1:
        sim_prefix = f"{sys.argv[1]}"
    else:
        sim_prefix = f"{date}"

    sim_name = "2D_eff_index_"
    sim_name += "cavity_" if cavity_parameters["N_rings"] > 0 else ""
    sim_name += "and_outcoupler_" if outcoupler_parameters["N_rings"] > 0 else ""
    sim_name += f"{sim_prefix}_Exy_"
    sim_name += f"angle{source_tilt*180/np.pi:.2f}_wv{1/f*1e3:.1f}"#"_n_eff_h{n_eff_h:.4f}"#point{Z_f:.0f}_


    sim = Simulation(sim_name,symmetries=[]) #mp.Mirror(mp.X),mp.Mirror(mp.Y)])# mp.Mirror(mp.Y,phase=-1) ])#mp.Mirror(mp.Y,phase=-1)])#
    sim.extra_space_xy += wavelength#/n_eff_l
    sim.eps_averaging = False
    sim.force_complex_fields = False
    sim.init_geometric_objects( eff_index_info = eff_index_info,
                                resolution = 40,
                                pattern_type = pattern_type,
                                cavity_parameters = cavity_parameters)

    if empty:
        sim.empty = True
        sim.name += '_empty'
    else:
        sim.empty = False

    sim.init_sources_and_monitors(f, df, source_pos=mp.Vector3(), #x=-sim.cavity_r_size - 0.1),
                                         source_tilt=source_tilt, allow_profile=False)# y=1e-3

    # raise Exception()1


    sim.init_sim()
    if mp.am_really_master():
        fig = plt.figure(dpi=150, figsize=(10,10))
        plot = sim.plot2D(eps_parameters={"interpolation":'none',"cmap":'gnuplot'})

        fig.colorbar(plot.images[0])
        fig.savefig(f'{sim.name}-xy.jpg')
        # plt.show()
        plt.close()

    # mp.verbosity(0)
    step_functions = []
    if sim.harminv_instance != None :
        step_functions.append( mp.after_time(150, sim.harminv_instance) )

    step_functions.append(mp.at_every(0.05,save_fields))
    sim.run(*step_functions, until=sim_end)
    if df == 0 :
        sim.run(save_fields, until=1/f * 5 ) # an integer number of periods

    print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to run\n')

    # plt.close()

    t = np.round(sim.round_time(), 2)

    data2save = {}
    if sim.harminv_instance != None :
        resonances_Q = []
        resonances_f = []
        for mode in  sim.harminv_instance.modes :
            if np.abs(mode.Q) > 10 :
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
        for j in range(len(sim.field_profile.freq)):
            data2save[f"field_profile_Hz_{j}"] = sim.get_dft_array(sim.field_profile, mp.Hz, j)
            data2save[f"field_profile_Ez_{j}"] = sim.get_dft_array(sim.field_profile, mp.Ez, j)
            data2save[f"field_profile_Ey_{j}"] = sim.get_dft_array(sim.field_profile, mp.Ey, j)
            data2save[f"field_profile_Ex_{j}"] = sim.get_dft_array(sim.field_profile, mp.Ex, j)
        data2save["field_profile_Eps"] = sim.get_array(mp.Dielectric,
                                                       center = sim.field_profile.regions[0].center,
                                                       size = sim.field_profile.regions[0].size)
        (x, y, _, _) = sim.get_array_metadata(center = sim.field_profile.regions[0].center,
                                              size = sim.field_profile.regions[0].size)
        data2save["field_profile_x"] = x
        data2save["field_profile_y"] = y
        data2save["field_profile_frequencies"] = np.array(sim.field_profile.freq)


    spectra = []
    for monitor in sim.spectrum_monitors :
        spectrum_f = np.array(mp.get_flux_freqs(monitor))
        spectra.append(np.array(mp.get_fluxes(monitor)))

    if len(sim.Ex) > 0:
        data2save["E_x"] = sim.Ex
        data2save["E_y"] = sim.Ey
        data2save["E_z"] = sim.Ez

    if len(spectra) > 0 :
        data2save["wavelength"] = 1/spectrum_f*1e3
        data2save["spectra"] = spectra
    if sim.field_FT != None :
        central_FT_x = np.array( [sim.get_dft_array(sim.field_FT, mp.Ex, j) for j in range(len(sim.field_FT.freq))] )
        central_FT_y = np.array( [sim.get_dft_array(sim.field_FT, mp.Ey, j) for j in range(len(sim.field_FT.freq))] )
        data2save["FT_x"] = central_FT_x
        data2save["FT_y"] = central_FT_y

    if len(data2save) > 0:
        mpo.savemat(f'{sim.name}_spectra_t{t}.mat', data2save)

    return data2save, sim.name


#%% geometry and simulation parameters
if __name__ == "__main__":              # good practise in parallel computing

    anisotropy = 0

    wavelength = .590

    period = .281 #round(wavelength/(n_eff_l+n_eff_h),3 )

    data = io.loadmat("Lumerical-Objects/multilayer_design/designs/TE_N7_dispersion_azoPPA_1.615.mat")
    n_eff = itp.RegularGridInterpolator((data["d"][0], data["lambda"][0]), data["n_eff"])
    # data = io.loadmat("16_29_50_WR_cavity_and_outcoupler_pos_280_D661_FF0d4_Ndbr30_Nout10_charge0_x1_15nm_1025C_RADIAL_compression.mat")

    # r = data["R"][0]
    # theta = data["theta"][0]
    # Z = data["Z"] + 60.8
    # theta = np.concatenate((theta, -theta[0:1]))
    # Z = np.concatenate((Z[r<9,:], Z[r<9,0:1]),axis=1)
    # r = r[r<9]
    # Z = Z[r>.1,:]
    # r = r[r>.1]
    # Z_interp = itp.RegularGridInterpolator((r, theta), Z)
    # R, THT = np.meshgrid(r,theta)
    # plt.figure()
    # ax = plt.subplot(111,projection='3d')
    # ax.plot_surface(R,THT, Z.transpose())

    # Z_f = lambda rr, tht: Z_interp((rr,tht))

    # data = io.loadmat("topo_resampled2.mat")
    # x = data["xx"][0]
    # y = data["yy"][0]
    # Z = data["topod"] + 65 - 4.2
    # Z_interp = itp.RegularGridInterpolator((y, x), Z)



    Z_f = lambda x, y: 1 # Z_interp((y,x))
    n_eff_h = n_eff([31e-9, wavelength*1e-6])[0]
    n_eff_l = n_eff([ 3e-9, wavelength*1e-6])[0]
    n_eff_h_v = [ n_eff_h ]#, 1.1045]
    n_eff_l_v = [ n_eff_l ]#, 1.0395]
    n_eff_mod_l = n_eff([15e-9, wavelength*1e-6])[0] - n_eff_l
    n_eff_mod_h = n_eff([39e-9, wavelength*1e-6])[0] - n_eff_h
    n_eff_spacer = n_eff([65e-9, wavelength*1e-6])[0]
    #%% load susceptibilities data.
    # Even though the variable are still called n_eff, they refer
    # to epsilon susceptibilities. mpo.Medium() can handle it

    # data = io.loadmat("bsw_lorentz_fit_idx1.000_th0.mat")
    # n_eff_l = [ a for a in data["optimal_fit_2"][0]]
    # data = io.loadmat("bsw_lorentz_fit_idx1.729_th23.mat")
    # n_eff_h = [ a for a in data["optimal_fit_2"][0]]

    #%%
    D = 0.560 #

    # crete input vector for parallell pool. It has to be a list of tuples,
    # where each element of the list represent one iteration and thus the
    # element of the tuple represent the inputs.
    empty = True
    tuple_list = [(wavelength,
                    n_eff_h_v[0], n_eff_l_v[0], n_eff_spacer,
                    D, period,
                    empty,
                    0,
                    n_eff_mod_l,
                    n_eff_mod_h )]

    empty = False

    j = 1
    j = 0           # resets  tiple list (insted of commenting all previous lines)
    tuple_list = []

    for source_tilt in np.linspace(np.pi/4, +np.pi/2, 1)[:]:

    # for source_pos in [0]: # 0, period/4, period/2]:

    # for i in range(len(n_eff_h_v)) :
    #     n_eff_h = n_eff_h_v[i]
    #     n_eff_l = n_eff_l_v[i]

    # for anisotropy in np.linspace(0,5, 1):

    # source_tilt = 0
    # for D in [1] : # np.linspace(0, 1, 50):


    # test the eff indeces extracted from the contour
    # of resonances as a function of thicknesses
    # data = io.loadmat("cross_points_7edc4aea0b.mat")
    # points584 = data["points"]
    # data = io.loadmat("cross_points_68fa746847.mat")
    # points596 = data["points"]

    # for i in range(1):#len(points584)):

    # test various thicknesses
    # for th_low in np.linspace(0, 65, 25):
    #     for th_high in np.linspace(0, 65, 25):

        for wavelength in np.linspace(.590, .615, 1):
            th = np.linspace(0,70,50)
            n_eff_tmp = itp.interp1d(th, n_eff( (th*1e-9, wavelength*1e-6*np.ones(50) ) ))
            n_eff_wv = lambda th : n_eff_tmp(th).item()
            n_eff_h      = n_eff_wv(32) # points584[i,1])
            n_eff_l      = n_eff_wv(3) # points584[i,0])
            n_eff_mod_l  = n_eff_wv(16) - n_eff_wv(3)# points596[i,0]) - n_eff_wv(points584[i,0])
            n_eff_mod_h  = n_eff_wv(41) - n_eff_wv(32)# points596[i,1]) - n_eff_wv(points584[i,1])
            n_eff_spacer = n_eff_wv(65)

            source_pos=0
            if n_eff_h > n_eff_l:
                tuple_list.append( (wavelength,
                                    n_eff_h, n_eff_l, n_eff_spacer,
                                    D, period,
                                    empty,
                                    source_pos, source_tilt,
                                    n_eff_mod_l,
                                    n_eff_mod_h, n_eff_wv, Z_f ) )
                j += 1

    mp.verbosity(1)
    # mp.quiet(True)
    output = []
    names = []
    tuple_list.reverse()
    t0 = time.time()
    try:
        from mpi4py import MPI
    except:
        non_parallel_conda = True
    else:
        non_parallel_conda = False

    if len(sys.argv) > 2:
        if sys.argv[2] == "parallel_grid":
            non_parallel_conda = True
        else:
            bash_parallel_run = (sys.argv[2] == "parallel_bash")

    if len(sys.argv) < 2 or non_parallel_conda :
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

    elif bash_parallel_run :
        N_jobs = int(sys.argv[-1])
        j = int(sys.argv[3])

        N_list = len(tuple_list)
        if N_list < N_jobs :
            raise ValueError(f"Number of jobs should be lower than number of loop iterations ({N_list})")

        reminder = np.mod(N_list,N_jobs)
        N_loops_per_job = int(N_list/N_jobs)
        if j < reminder:
            N_loops_per_job += 1

        data_list = []
        name_list = []
        for i in range(N_loops_per_job):
            t1 = time.time()
            if j < reminder:
                tuple_index = j*N_loops_per_job + i
            else:
                tuple_index = reminder*(N_loops_per_job+1) + (j-reminder)*N_loops_per_job + i

            if tuple_index >= N_list :
                continue
            data, name = run_parallel(*tuple_list[tuple_index])
            # data_list.append(data)
            # name_list.append(name)
            print(f'It has run for {convert_seconds(time.time()-t1)}, {i+1}/{N_loops_per_job}')
            print(f'It will take roughly {convert_seconds((time.time()-t0)/(i+1)*(N_loops_per_job-i-1))} more')

    else:
        # mp.reset_meep()
        comm = MPI.COMM_WORLD
        N_jobs = int(sys.argv[-1])
        print(f'number of processor is {mp.count_processors()}')
        j = mp.divide_parallel_processes(N_jobs)

        N_list = len(tuple_list)
        if N_list < N_jobs :
            raise ValueError(f"Number of jobs should be lower than number of loop iterations ({N_list})")

        reminder = np.mod(N_list,N_jobs)
        N_loops_per_job = int(N_list/N_jobs)
        if j < reminder:
            N_loops_per_job += 1

        data_list = []
        name_list = []
        for i in range(N_loops_per_job):
            t1 = time.time()
            if j < reminder:
                tuple_index = j*N_loops_per_job + i
            else:
                tuple_index = reminder*(N_loops_per_job+1) + (j-reminder)*N_loops_per_job + i

            if tuple_index >= N_list :
                continue
            data, name = run_parallel(*tuple_list[tuple_index])
            # data_list.append(data)
            # name_list.append(name)
            print(f'It has run for {convert_seconds(time.time()-t1)}, {i+1}/{N_loops_per_job}')
            print(f'It will take roughly {convert_seconds((time.time()-t0)/(i+1)*(N_loops_per_job-i-1))} more')

        # if mp.am_really_master():
        #     output.extend(data_list)
        #     names.extend(name_list)
        #     for src in range(1, N_jobs):
        #         output.extend( comm.recv(source=src, tag=11) )
        #         names.extend ( comm.recv(source=src, tag=12) )
        #         # comm.recv(source=src, tag=11)
        #         # comm.recv(source=src, tag=12)
        # else:
        #     comm.send(data_list, dest=0, tag=11)
        #     comm.send(name_list, dest=0, tag=12)
        #     exit()
    print(f'Total took {convert_seconds(time.time()-t0)}')

    #%%
    # plt.figure()
    # wv = output[0]["wavelength"]
    # s0 = output[0]["spectra"][0]
    # s1 = output[1]["spectra"][0]/s0
    # s2 = output[2]["spectra"][0]/s0
    # s3 = output[3]["spectra"][0]/s0
    # plt.semilogy(wv, s1, wv, s2, wv, s3)
    # plt.grid(True)
    # plt.xlabel("wavelength")
