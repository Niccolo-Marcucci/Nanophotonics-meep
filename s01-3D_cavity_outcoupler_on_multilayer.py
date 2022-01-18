#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import meep as mp
import numpy as np
# import h5py as h5
#import scipy as sp
from scipy import optimize as op
from scipy import interpolate as itp
from matplotlib import pyplot as plt
# from multiprocessing import Pool
# from mpl_toolkits.mplot3d import Axes3D
import meep_objects as mpo
# import io
import sys
import json
import time


# from mayavi import mlab

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

    def __init__(self, sim_name='simulation_buried', dimensions=3, symmetries = []):

        self.name = sim_name

        self.extra_space_xy = .3

        self.PML_width = .3

        self.z_top_air_gap = 0.7

        self.substrate_thickness = .2

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
        try:
            if self._empty :
                self.geometry = self._empty_geometry
            else:
                self.geometry = self._empty_geometry
                self.geometry.extend(self._geometry)
        except AttributeError:
            raise AttributeError("cannot assign 'empty' property before initializing the geometry")
        self.reset_meep()

    def init_geometric_objects(self, multilayer_file, used_layer=-3, res_scaling=1, use_BB=True,
                               pattern_type='positive', cavity_parameters={}, outcoupler_parameters={}):

                               # D=5, grating_period=0.2, N_rings=10, N_arms=0, lambda_bsw=0.4,
                               # scatter_length=0.4, scatter_width=0.1, scatter_tilt=0,
                               # scatter_shape='', scatter_disposition='filled', topology='spiral', pattern_type='positive') :
        self._geometry = []
        self._empty_geometry = []

        self.cavity_r_size = (cavity_parameters["D"]/2 + cavity_parameters["period"] * cavity_parameters["N_rings"]) * (cavity_parameters["N_rings"]>0)
        self.outcou_r_size = (outcoupler_parameters["D"]/2 + outcoupler_parameters["period"] * outcoupler_parameters["N_rings"]) * (outcoupler_parameters["N_rings"]>0)

        self.domain_x = self.domain_y = 2*(self.cavity_r_size + self.outcou_r_size + self.extra_space_xy)

        multilayer, multilayer_thickness, design_specs = mpo.dielectric_multilayer(
            design_file = multilayer_file,
            substrate_thickness = self.substrate_thickness + .5 + self.PML_width,
            x_width = self.domain_x + .5 + 2*self.PML_width,
            y_width = self.domain_y + .5 + 2*self.PML_width,
            unit = 'um',
            exclude_last_layer = False,
            buried = True)

        print(design_specs)
        self._empty_geometry.extend(multilayer)            # keep multilayer even if empty

        if pattern_type == 'positive':
            grating_index = np.real(design_specs['idx_layers'][used_layer])

            dummy_layer = mp.Block(
                material = mp.Medium(index = np.real(design_specs['idx_layers'][used_layer+1])),
                size     = mp.Vector3(self.domain_x + .5 + 2*self.PML_width,
                                      self.domain_y + .5 + 2*self.PML_width,
                                      design_specs['d_layers'][used_layer]),
                center   = mp.Vector3(0, 0, 0))
            self._empty_geometry.append(dummy_layer)       # part of the multilayer

        elif pattern_type == 'negative':
            grating_index = np.real(design_specs['idx_layers'][used_layer+1])

        else :
            raise ValueError(f'patter type "{pattern_type}" is unknown')

        if cavity_parameters["N_rings"] > 0:
            cavity = mpo.spiral_grating(
                medium_groove = mp.Medium(index=grating_index),
                D = cavity_parameters["D"],
                FF = cavity_parameters["FF"],
                DBR_period = cavity_parameters["period"],
                N_rings = cavity_parameters["N_rings"],
                N_arms = 0,
                thickness = float(design_specs['d_layers'][used_layer]),
                center = mp.Vector3())
            self._geometry.extend(cavity)

        if outcoupler_parameters["type"] == 'pol_splitting' and outcoupler_parameters["N_rings"] > 0:
            outcoupler = mpo.pol_splitting_grating(
                medium_groove = mp.Medium(index=grating_index),
                D = self.cavity_r_size*2 + outcoupler_parameters["D"],
                metasurface_period = outcoupler_parameters["period"],
                scatter_length = outcoupler_parameters["scatter_length"],
                scatter_width = outcoupler_parameters["scatter_width"],
                scatter_tilt = outcoupler_parameters["scatter_tilt"],
                scatter_shape = outcoupler_parameters["scatter_shape"],
                scatter_disposition = outcoupler_parameters["scatter_disposition"],
                topology = outcoupler_parameters["topology"],
                n_rings = outcoupler_parameters["N_rings"],
                n_arms = outcoupler_parameters["N_arms"],
                lambda_bsw = outcoupler_parameters["lambda_bsw"],
                thickness = float(design_specs['d_layers'][used_layer]),
                center = mp.Vector3())
            self._geometry.extend(outcoupler)

        elif outcoupler_parameters["type"] == 'spiral' and outcoupler_parameters["N_rings"] > 0:
            outcoupler = mpo.spiral_grating(
                medium_groove = mp.Medium(index=grating_index),
                D = self.cavity_r_size*2 + outcoupler_parameters["D"],
                FF = outcoupler_parameters["FF"],
                DBR_period = outcoupler_parameters["period"],
                N_rings = outcoupler_parameters["N_rings"],
                N_arms = outcoupler_parameters["N_arms"],
                thickness = float(design_specs['d_layers'][used_layer]),
                center = mp.Vector3())
            self._geometry.extend(outcoupler)

        if use_BB:
            beam_block = mp.Cylinder(
                            center = mp.Vector3(0, 0, self.z_top_air_gap-0.1),
                            radius = 3/4 * (self.cavity_r_size + outcoupler_parameters["D"]/2),
                            height = 0.02,
                            material = mp.metal)
            beam_block.name = 'Beam_block'

            self._geometry.append(beam_block)

        self.empty = False                  # this  will add all geometric objects to the simulation

        self.domain_z = self.substrate_thickness + multilayer_thickness + self.z_top_air_gap

        # resolution is 10 points per wavelength in the highest index material time a scale factor
        self.resolution = int(10/(1/f/np.real(np.max(design_specs['idx_layers']))) * res_scaling)

        self.name = self.name + f'res{sim.resolution}'
        self.filename_prefix = self.name

        # round domain with an integer number of grid points
        self.grid_step = 1/self.resolution

        self.cell_size = mp.Vector3(self.domain_x + 2*self.PML_width,
                                    self.domain_y + 2*self.PML_width,
                                    self.domain_z + 2*self.PML_width)
        # make domain an integer number of voxels
        Nx = int(self.cell_size.x / self.grid_step)
        Nx -= np.mod(Nx,2) + 1      # make odd
        self.cell_size.x = Nx * self.grid_step
        Ny = int(self.cell_size.y / self.grid_step)
        Ny -= np.mod(Ny,2) + 1
        self.cell_size.y = Ny * self.grid_step
        Nz = int(self.cell_size.z / self.grid_step)
        Nz -= np.mod(Nz,2) + 1
        self.cell_size.z = Nz * self.grid_step

        print()
        print(f"Number of voxels is ({Nx}x{Ny}x{Nz}) = {Nx*Ny*Nz/1e6} Mln")
        print(f"Minimum expected memory is {96*Nx*Ny*Nz/2**30:.2f}GB")
        print()

        self.geometry_center = mp.Vector3(0, 0, -(self.cell_size.z/2 - self.z_top_air_gap - self.PML_width - np.sum(design_specs['d_layers'][used_layer+1:-1]) - design_specs['d_layers'][used_layer]/2))

        self.boundary_layers = [mp.PML(self.PML_width)]
        # print( [self.cell_size.x / self.

        with open(f'{self.name}.json', 'w') as fp:
            data2save = {"multilayer": multilayer_file,
                         "pattern_type": pattern_type,
                         "use_beam_block": use_BB}

            if cavity_parameters["N_rings"] > 0:
                data2save["cavity_parameters"] = cavity_parameters

            if outcoupler_parameters["N_rings"] > 0:
                data2save["outcoupler_parameters"] = outcoupler_parameters

            json.dump(data2save, fp,  indent=4)


    def init_sources_and_monitors(self, f, df, allow_farfield=True) :
        self.sources = [ mp.Source(
            src = mp.ContinuousSource(f,fwidth=0.1) if df==0 else mp.GaussianSource(f,fwidth=df),
            center = mp.Vector3(),
            size = mp.Vector3(),
            component = mp.Ez)]

        self.nearfield_monitor = None
        self.harminv_instance = None
        self.spectrum_monitors = []

        if self.outcou_r_size > 0 and allow_farfield :
            nearfield = mp.Near2FarRegion(
                center = mp.Vector3(0, 0, self.z_top_air_gap - 0.03),
                size = mp.Vector3(self.domain_x-2*self.grid_step, self.domain_y-2*self.grid_step, 0),
                direction = mp.Z)

            self.nearfield_monitor = self.add_near2far(f, 0, 1, nearfield)#, yee_grid=True))

        if self.cavity_r_size > 0 :
            nfreq = 1000
            fluxr = mp.FluxRegion(
                center = mp.Vector3(0, self.cavity_r_size, 0),
                size = mp.Vector3(.1,0,0),
                direction = mp.Y)
#            self.spectrum_monitors.append(self.add_flux(f, df, nfreq, fluxr))#, yee_grid=True))

            if not self.empty:
                self.harminv_instance =  mp.Harminv(mp.Hx, mp.Vector3(), f, df)




#%% geometry and simulation parameters
c0 = 1
wavelength = 0.570
wwidth = .20
f = c0 / wavelength

fmax = c0 / (wavelength - wwidth/2)
fmin = c0 / (wavelength + wwidth/2)
df = fmax - fmin

sim_end = 1e-4

n_eff_l = 1.6642
n_eff_h = 1.7899
n_eff_FF0d5 = n_eff_h*.5 + n_eff_l*.5


file = 'design_TM_gd3_buriedDBR_onSiO2'
buried = True
pattern_type = 'positive'           # 'positive' or 'negative'
out_grating_type = 'spiral'         # 'spiral' or 'polSplitting' or 'only'

# cavity info
N_cavity = 30
cavity_period = .165 # wavelength / n_eff_FF0d5 / 2
D_cavity = .400 # cavity_period * 1.4

# pol splitting info
FF_pol_splitter = .3
FF = FF_pol_splitter
n_eff = n_eff_h*(1-FF) + n_eff_l*FF if pattern_type=='positive' else n_eff_h*FF + n_eff_l*(1-FF)
scatter_disposition='filled'        # 'radial' or 'filled'
D_phi = np.pi/3;
sigma = -1;                         # select for circl left or circ right
K_bsw = 2*np.pi * n_eff / wavelength
m = 1                               # ordinary grating order
s = (m*2*np.pi + sigma * 2*D_phi) / K_bsw
outcoupler_period = s

# outcoupler info
N_outcoupler = 0
d_cavity_out = .5
charge = 3

cavity_parameters = {
    "D": D_cavity,
    "FF": .5,
    "period": cavity_period,
    "N_rings": N_cavity}

spiral_parameters = {
    "type": 'spiral',
    "D": d_cavity_out,
    "FF": .5,
    "period": wavelength / n_eff_FF0d5,
    "N_rings": N_outcoupler if out_grating_type=='spiral' else 0,
    "N_arms": charge}

polSplitter_parameters = {
    "type": 'pol_splitting',
    "D": d_cavity_out,
    "period": outcoupler_period,
    "scatter_length": outcoupler_period*0.9,
    "scatter_width": outcoupler_period*FF_pol_splitter,
    "scatter_tilt": D_phi,
    "scatter_shape": '',
    "scatter_disposition": scatter_disposition,
    "topology": 'spiral',
    "N_rings": N_outcoupler if out_grating_type=='polSplitting' else 0,
    "N_arms": charge,
    "lambda_bsw": wavelength/n_eff,
    "sigma": sigma}

t0 = time.time()


date = time.strftime('%y%m%d-%H%M%S')#'211001-121139'#
if len(sys.argv) > 1:
    sim_prefix = f"{sys.argv[1]}"
else:
    sim_prefix = f"{date}"

sim_name = "cavity_" if N_cavity > 0 else ""
sim_name += f"{out_grating_type}_" if N_outcoupler > 0 else ""
sim_name += f"{sim_prefix}_{file}"
sim_name += f"_charge{charge}" if N_outcoupler > 0 else ""

# sim_name += f"_{parameter_to_loop}"

sim = Simulation(sim_name)
sim.extra_space_xy += wavelength/n_eff_l * (charge > 0)
sim.eps_averaging = True

sim.init_geometric_objects( multilayer_file = f"./Lumerical-Objects/multilayer_design/designs/{file}",
                            used_layer = -3 if buried else -2,
                            res_scaling = 1,
                            use_BB = False,
                            pattern_type = pattern_type,
                            cavity_parameters = cavity_parameters,
                            outcoupler_parameters = spiral_parameters if out_grating_type=='spiral' else polSplitter_parameters)
sim_name = sim.name
sim.init_sources_and_monitors(f, df)
mp.verbosity(2)
# mpo.create_openscad(sim,1000)
sim.init_sim()

# raise ValueError()

print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to initiate\n')

#%%
simsize = sim.cell_size
center  = sim.geometry_center

max_epsilon = 2.53**2

fig = plt.figure(dpi=200)
plot = sim.plot2D( output_plane=mp.Volume(center=center, size=mp.Vector3(0,simsize.y,simsize.z)),
                   labels=True,
                   eps_parameters={"interpolation":'none',"cmap":'gnuplot', "vmin":'0.5', "vmax":max_epsilon} )
try:
    fig.colorbar(plot.images[0], orientation="horizontal")
except:
    plt.close()
    print("Only one of the parallel jobs jobs will print the image")
else:
    fig.savefig(f'{sim_name}_section-yz.jpg')
    plt.close()

fig = plt.figure(dpi=200)
plot = sim.plot2D( output_plane=mp.Volume(center=mp.Vector3(z=-.00), size=mp.Vector3(simsize.x,simsize.y)),
                   labels=True,
                   eps_parameters={"interpolation":'none',"cmap":'gnuplot', "vmin":'0.5', "vmax":max_epsilon})
try:
    fig.colorbar(plot.images[0])
except:
    plt.close()
    print("Only one of the parallel jobs jobs will print the image")
else:
    fig.savefig(f'{sim_name}_section-xy.jpg')
    plt.close()

# sim.output_epsilon(f'{sim_name}_eps')
# eps_data = sim.get_epsilon()
# mpo.savemat(f'{sim_name}_eps.mat', {"eps_data": eps_data})
# x, y, z, w = [np.array(tmp) for tmp in sim.get_array_metadata()]
# mpo.plot_image(z, y, eps_data[:,:,84], vmax=9.0, vmin=1.0)
# mpo.plot_image(y, z, eps_data[int(eps_data.shape[0]/2)+1,:,:])#, vmax=9.0, vmin=-1.0)

# mpo.plot_data_section(eps_data)

# # s = mlab.con(x,y,z,w)=sim.get_array_metadata()tour3d(eps_data, colormap="YlGnBu")
# # mlab.show()

#%%
# raise RuntimeError("comment this line to run til the end")
def print_time(sim):
    print(f'\n\nSimulation is at {sim.round_time()} \n It has run for {convert_seconds(time.time()-t0)}\n')

t0 = time.time()
mp.verbosity(1)
for i in range(1):

    # fig = plt.figure(dpi=100)
    # Animate = mp.Animate2D( sim, fields=mp.Ez, f=fig, realtime=False, normalize=True,
    #                         output_plane=mp.Volume(center=mp.Vector3(), size=mp.Vector3(simsize.x,simsize.y,0)),
    #                         eps_parameters={"interpolation":'none',"vmin":'0'})

    # sim.run(mp.at_every(.1, Animate),until=30)
    # Animate.to_mp4(10,f'{sim_name}_section.mp4')

    step_functions = [mp.at_every(5,print_time)]
    if sim.harminv_instance != None :
        step_functions.append( mp.after_sources(sim.harminv_instance) )


    sim.run(*step_functions, until=100)# until_after_sources=mp.stop_when_fields_decayed(1, mp.Ez, mp.Vector3(), 1e-1))
    # sim.run(until_after_sources=mp.stop_when_dft_decayed(minimum_run_time=10))

    # Animate.to_mp4(10,f'{sim_name}_section.mp4')
    print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to run\n')

    t = np.round(sim.round_time(), 2)

    if sim.nearfield_monitor != None :
        ex_near, ey_near = [sim.get_dft_array(sim.nearfield_monitor, field, 0) for field in [mp.Ex, mp.Ey]]
        mpo.savemat(f'{sim_name}_nearfield_t{t}.mat', {'Ex': ex_near, 'Ey': ey_near,
                                                       'Lx': sim.nearfield_monitor.regions[0].size.x,
                                                       'Ly': sim.nearfield_monitor.regions[0].size.y})
    if sim.harminv_instance != None :
        resonances_Q = []
        resonances_f = []
        for mode in  sim.harminv_instance.modes :
            if np.abs(mode.Q) > 0 :
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
            resonance_table.append([np.round(1/resonances_f[l]*1e3, 1), np.int(resonances_Q[l])] )
        if N_resonances == 0 :
            resonance_table.append([ 0, 0 ])
        print()
        print(resonance_table)
        print()

        with open(f'{sim_name}_output.json', 'a') as fp:
            data2save = {f"resonance_table_t{t}": resonance_table}
            json.dump(data2save, fp,  indent=4)

spectra = []
for monitor in sim.spectrum_monitors :
    spectrum_f = np.array(mp.get_flux_freqs(monitor))
    spectra.append(np.array(mp.get_fluxes(monitor)))

if len(spectra) > 0 :
    sim.empty = True
    sim.init_sources_and_monitors(f, df, allow_farfield=False)

    sim.run(mp.at_every(5,print_time), until=t)

    spectra_out = []
    for i, monitor in enumerate(sim.spectrum_monitors) :
        spectrum_empty = mp.get_fluxes(monitor)
        spectra_out.append( np.array(spectrum_empty) )

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111)
    plt.plot(1/spectrum_f*1e3, spectra_out[0])
    # plt.xlim(540,660)
    # plt.ylim(-2,2)
    ax.grid(True)
    plt.xlabel('wavelength')
    plt.ylabel('Transmission')
    ax2 = fig.add_subplot(336)
    # plt.title('Table of the resonances')
    collabel=[ "Wavelength", "Quality"]
    rowlabel=[ f'{i}' for i in range(len(resonance_table))]
    ax2.axis('tight')
    ax2.axis('off')
    the_table = ax2.table(cellText=resonance_table, colLabels=collabel, rowLabels=rowlabel,loc='center')

    fig.savefig(f'{sim_name}_spectrum_cavity.jpg')
    plt.close(fig)

    mpo.savemat(f'{sim_name}_spectra_t{t}.mat', {"wavelength": 1/spectrum_f*1e3,
                                                  "spectra"   : spectra,
						  "spectrum_empty": spectra_out,
                                                  "resonances" : resonance_table})
