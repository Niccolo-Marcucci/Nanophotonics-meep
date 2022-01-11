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
# import sys
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



class Simulation():

    def __init__(self, sim_name='simulation_buried', dimensions=2, symmetries = [], empty = False):

        self.name = sim_name

        self.extra_space_xy = .3

        self.PML_width = .5

        self.z_top_air_gap = 0.7

        self.substrate_thickness = .2

        self._empty = empty

        self.sim = mp.Simulation(
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

        if self._empty :
            self.sim.geometry = []
        else:
            self.sim.geometry = self.geometry

        self.sim.reset_meep()

    def init_geometric_objects(self, multilayer_file, D=5, grating_period=0.2, N_rings=10, N_arms=0, lambda_bsw=0.4,
                               pattern_type = 'positive') :
        self.geometry = []

        self.domain_x = grating_period*N_rings*2 + D + self.extra_space_xy*2

        multilayer, multilayer_thickness, design_specs = mpo.dielectric_multilayer(
                        design_file = multilayer_file,
                        substrate_thickness = self.substrate_thickness + .5 + 2*self.PML_width,
                        x_width = self.domain_x + .5 + 2*self.PML_width,
                        unit = 'um',
                        exclude_last_layer = False,
                        buried = True,
                        axis = "y")
        print(design_specs)

        self.geometry.extend(multilayer)

        beam_block = mp.Cylinder(
                        center = mp.Vector3(0, self.z_top_air_gap-0.1,0),
                        radius = D/2,
                        height = 0.02,
                        material = mp.metal,
                        axis = mp.Vector3(0, 1, 0))
        beam_block.name = 'Beam_block'

        self.geometry.append(beam_block)

        # if pattern_type == 'positive':
        #     grating_index = np.real(design_specs['idx_layers'][-3])

        #     dummy_layer = mp.Block(material = mp.Medium(index = np.real(design_specs['idx_layers'][-2])),
        #                             size     = mp.Vector3(self.domain_x+.5 + 2*self.PML_width,
        #                                                   design_specs['d_layers'][-3],1),
        #                             center   = mp.Vector3(0, 0, 0))
        #     self.geometry.append(dummy_layer)

        # elif pattern_type == 'negative':
        #     grating_index = np.real(design_specs['idx_layers'][-2])

        # else :
        #     raise ValueError(f'patter type "{pattern_type}" is unknown')

        # outcoupler = mpo.spiral_grating(
        #                 medium_groove = mp.Medium(index=grating_index),
        #                 D = D,
        #                 FF = 0.5,
        #                 DBR_period = grating_period,
        #                 N_rings = N_rings,
        #                 N_arms = N_arms,
        #                 thickness = design_specs['d_layers'][-3],
        #                 orientation = mp.Vector3(0, 1, 0))

        # self.geometry.extend(outcoupler)

        if not self.empty:
            self.sim.geometry = self.geometry

        self.domain_y = self.substrate_thickness + multilayer_thickness + self.z_top_air_gap

        # resolution is 10 points per wavelength in the highest index material time a scale factor
        self.sim.resolution = int(10/(1/f/np.real(np.max(design_specs['idx_layers']))) * 4)
        print(1/self.sim.resolution)

        # round domain with an integer number of grid points
        self.grid_step = 1/self.sim.resolution
        # self.domain_x = int(self.domain_x/self.grid_step) * self.grid_step
        # self.domain_y = int(self.domain_y/self.grid_step) * self.grid_step
        # self.PML_width = int(self.PML_width/self.grid_step) * self.grid_step
        # print(self.PML_width)
        # print(self.grid_step)
        # self.domain_z = int(self.domain_z/self.grid_step) * self.grid_step

        self.sim.cell_size = mp.Vector3(self.domain_x + 2*self.PML_width,
                                        self.domain_y + 2*self.PML_width)

        self.sim.geometry_center = mp.Vector3(0, -(self.sim.cell_size.y/2 - self.z_top_air_gap - self.PML_width - design_specs['d_layers'][-2] - design_specs['d_layers'][-3]/2), 0)

        self.sim.boundary_layers = [mp.PML(self.PML_width)]
        # print( [self.sim.cell_size.x / self.sim.

    def init_sources_and_monitors(self, f, df, k) :

        self.sim.sources = [
                            mp.Source(src=mp.ContinuousSource(f,fwidth=0.1) if df==0 else mp.GaussianSource(f,fwidth=df),
                                  center = mp.Vector3(),
                                  size = mp.Vector3(),
                                  component=mp.Ey)]
                            # mp.EigenModeSource(src=mp.ContinuousSource(f,fwidth=0.1) if df==0 else mp.GaussianSource(f,fwidth=df),
                            #       center = mp.Vector3(-self.sim.cell_size.x/4)+self.sim.geometry_center,
                            #       size = mp.Vector3(y=self.sim.cell_size.y-2*self.PML_width),
                            #       direction=mp.AUTOMATIC,
                            #       eig_kpoint=k,
                            #       eig_band=mode_N,
                            #       eig_parity=mp.NO_PARITY,
                            #       eig_match_freq=True)]
                            # *mpo.plane_wave_source(f, df, k, center=mp.Vector3(0,.5,0),
                            #                    size=mp.Vector3(self.sim.cell_size.x),
                            #                    inc_plane_norm=mp.Vector3(0,0,1))]


        monitor_distance  = self.z_top_air_gap - 0.03

        nfreq = 1000

        self.monitors = []

        nearfield = mp.Near2FarRegion(mp.Vector3(0, monitor_distance, 0),
                                      size = mp.Vector3(self.domain_x-2*self.grid_step,1, 0),
                                      direction = mp.Y)

        # flux = self.sim.add_flux(f, 0, 1, mp.FluxRegion(center = mp.Vector3(+self.sim.cell_size.x/4)+self.sim.geometry_center,
        #                                                size = mp.Vector3(0,self.sim.cell_size.y-2*self.PML_width)))

        self.monitors.append(self.sim.add_near2far(f, 0, 1, nearfield))#, yee_grid=True))


    def create_openscad(self, scale_factor=1):
        try:
            import openpyscad as ops
        except ModuleNotFoundError:
            print("WARNING openpyscad is not installed in this environment")
        else:
            scad_name = f"{self.name}.scad"
            with open(scad_name,"w") as file:
                file.write(f"//Simulation {self.name}\n")

            for obj in self.geometry:
                if obj.__class__ == mp.Block:
                    cube = ops.Cube([obj.size.x*scale_factor,
                                      obj.size.y*scale_factor,
                                      obj.size.z*scale_factor], center = True)
                    tilt = np.arctan2(obj.e1.y, obj.e1.x)
                    if tilt != 0:
                        cube = cube.rotate([0, 0, tilt/np.pi*180])

                    index = np.round(np.sqrt(obj.material.epsilon_diag.x),2)
                    if index == 2.53:
                        color = [0, 0, 1]
                    elif index == 1.65:
                        color = [0, 1, 0]
                    elif index == 1.46:
                        color = [1, 0, 0]
                    elif index == 1.48:
                        color = [1, 1, 0]
                    elif index == 2.08:
                        color = [0, 1, 1]
                    cube = cube.color(color)

                    scad = cube

                elif obj.__class__ == mp.Cylinder:
                    cyl = ops.Cylinder(h=obj.height*scale_factor,
                                        r=obj.radius*scale_factor, center = True)

                    scad = cyl

                scad = scad.translate([obj.center.x*scale_factor,
                                       obj.center.y*scale_factor,
                                       obj.center.z*scale_factor])
                scad.write(scad_name, mode='a')

            sim_domain = ops.Cube([(self.sim.cell_size.x - self.PML_width) * scale_factor,
                                    (self.sim.cell_size.y - self.PML_width) * scale_factor,
                                    (self.sim.cell_size.z - self.PML_width) * scale_factor], center = True)
            sim_domain = sim_domain.color([.5, .5, .5, .5])
            sim_domain = sim_domain.translate([self.sim.geometry_center.x*scale_factor,
                                                self.sim.geometry_center.y*scale_factor,
                                                self.sim.geometry_center.z*scale_factor])


    # define aliases
    def run(self,*args,**kwargs):
        self.sim.run(*args,**kwargs)

    def init_sim(self):
        self.sim.init_sim()

    def round_time(self):
        self.sim.round_time()


#%% geometry and simulation parameters
c0 = 1
wavelength = 0.570
wwidth = .25
f = c0 / wavelength

fmax = c0 / (wavelength - wwidth/2)
fmin = c0 / (wavelength + wwidth/2)
df =  fmax - fmin

sim_end = 1e-4

n_eff_l = 1.6642
n_eff_h = 1.7899
n_eff = n_eff_h*.5 + n_eff_l*.5

pattern_type = 'positive'

outcoupler_period = round(wavelength/n_eff,3)
N_periods = 0
D = 5
charge = 0

t0 = time.time()

k0 = 2*np.pi / wavelength

kx = k0 * n_eff

ky = 0 #np.sqrt(k0**2 - kx**2, dtype=complex)

k = mp.Vector3(kx, ky, 0)
mode_N = 12
# file = 'design_TE_StellaEtAll_2019'
file = 'design_TM_gd3_buriedDBR_onSiO2'

sim_name = f"coupling_test_{file}_TMdipole"
sim = Simulation(sim_name)
# sim.empty = True
sim.init_geometric_objects(
                multilayer_file = f"Lumerical-Objects/multilayer_design/designs/{file}",
                D = D,
                grating_period = outcoupler_period,
                N_rings = N_periods,
                N_arms  = charge,
                pattern_type=pattern_type)

sim.init_sources_and_monitors(f, df, k)
mp.verbosity(1)
sim.init_sim()

# sim.create_openscad(scale_factor = 1e3)
# raise ValueError()

date = time.strftime('%y%m%d-%H%M%S')#'211001-121139'#
sim_suffix = f'res{sim.sim.resolution}_{date}'

print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to initiate\n')
#%%
simsize = sim.sim.cell_size
center  = sim.sim.geometry_center
fig = plt.figure(dpi=300)
sim.sim.plot2D( output_plane=mp.Volume(size=mp.Vector3(simsize.x,2*simsize.y)),
                labels=True,
                eps_parameters={"interpolation":'none',"cmap":'gnuplot', "vmin":'0'})
# fig.savefig(f'{sim_name}-{sim_suffix}_section.jpg')
# plt.close()
plt.show()


# plt.show(block=False)
# sim.sim.output_epsilon(f'{sim_name}_eps')
# eps_data = sim.sim.get_epsilon()
# mpo.savemat(f'{sim_name}_eps.mat', {"eps_data": eps_data})
# x, y, z, w = [np.array(tmp) for tmp in sim.sim.get_array_metadata()]
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
#
# sim.run(until=7)
# times_list = []
# for i in range(1):
#     sim.run(until=10)
#     # sim.run(until_after_sources=mp.stop_when_fields_decayed(1, mp.Ez, mp.Vector3(), sim_end))
#     t = np.round(sim.sim.round_time(), 5)
#     # sim.sim.save_near2far(near2far=sim.monitors[0], fname=f'{sim_suffix}_nearfield_t{t}')
#     times_list.append(t)

#     # ex_near, ey_near = [sim.sim.get_dft_array(sim.monitors[0], field, 0) for field in [mp.Ex, mp.Ey]]
#     # mpo.savemat(f'{sim_name}-{sim_suffix}_nearfield_t{t}.mat', {'Ex': ex_near, 'Ey': ey_near})

#     print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to run\n')

# flux = sim.monitors[0]
# res = sim.sim.get_eigenmode_coefficients(flux,[1], eig_parity=mp.NO_PARITY)
# incident_coeffs = res.alpha
# incident_flux = mp.get_fluxes(flux)
# incident_flux_data = sim.sim.get_flux_data(flux)

# Prepare the animator and record the steady state response
f = plt.figure(dpi=100)
# Animate = mp.Animate2D(sim.sim, fields=mp.Hz, f=f, realtime=False, normalize=False,
#                 eps_parameters={"interpolation":'none',"vmin":'0'})
sim.run(mp.at_every(0.1,print_time),until=2)


#%% Process the animation and view it
# filename = f'{sim_name}-{sim_suffix}_section.mp4'
# Animate.to_mp4(10,filename)

print('\a')
