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

    def __init__(self, sim_name='simulation_buried', dimensions=3, symmetries = [], empty = False):

        self.name = sim_name

        self.extra_space_xy = .3

        self.PML_width = .3

        self.z_top_air_gap = 0.7

        self.substrate_thickness = .2

        self._empty = empty

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

        if self._empty :
            self.geometry = []
        else:
            self.geometry = self._geometry

        self.reset_meep()

    def init_geometric_objects(self, multilayer_file, D=5, grating_period=0.2, N_rings=10, N_arms=0, lambda_bsw=0.4,
                               pattern_type = 'positive') :
        self._geometry = []

        self.domain_x = self.domain_y = grating_period*N_rings*2 + D + self.extra_space_xy*2

        multilayer, multilayer_thickness, design_specs = mpo.dielectric_multilayer(
                        design_file = multilayer_file,
                        substrate_thickness = self.substrate_thickness + .5 + 2*self.PML_width,
                        x_width = self.domain_x + .5 + 2*self.PML_width,
                        y_width = self.domain_y + .5 + 2*self.PML_width,
                        unit = 'um',
                        exclude_last_layer = False,
                        buried = True)
        print(design_specs)
        self._geometry.extend(multilayer)

        if pattern_type == 'positive':
            grating_index = np.real(design_specs['idx_layers'][-3])

            dummy_layer = mp.Block(material = mp.Medium(index = np.real(design_specs['idx_layers'][-2])),
                                    size     = mp.Vector3(self.domain_x+.5 + 2*self.PML_width,
                                                          self.domain_y+.5 + 2*self.PML_width,
                                                          design_specs['d_layers'][-3]),
                                    center   = mp.Vector3(0, 0, 0))
            self._geometry.append(dummy_layer)

        elif pattern_type == 'negative':
            grating_index = np.real(design_specs['idx_layers'][-2])

        else :
            raise ValueError(f'patter type "{pattern_type}" is unknown')


        outcoupler = mpo.spiral_grating(
                        medium_groove = mp.Medium(index=grating_index),
                        D = D,
                        FF = 0.5,
                        DBR_period = grating_period,
                        N_rings = N_rings,
                        N_arms = N_arms,
                        thickness = design_specs['d_layers'][-3],
                        orientation = mp.Vector3(0, 0, 1))

        self._geometry.extend(outcoupler)

        beam_block = mp.Cylinder(
                        center = mp.Vector3(0, 0, self.z_top_air_gap-0.1),
                        radius = D/2*3/4,
                        height = 0.02,
                        material = mp.metal)
        beam_block.name = 'Beam_block'

        self._geometry.append(beam_block)

        if not self.empty:
            self.geometry = self._geometry

        self.domain_z = self.substrate_thickness + multilayer_thickness + self.z_top_air_gap

                           # resolution is 10 points per wavelength in the highest index material time a scale factor
        self.resolution = int(10/(1/f/np.real(np.max(design_specs['idx_layers']))) * 1.5)
        print(1/self.resolution)

        # round domain with an integer number of grid points
        self.grid_step = 1/self.resolution
        # self.domain_x = int(self.domain_x/self.grid_step) * self.grid_step
        # self.domain_y = int(self.domain_y/self.grid_step) * self.grid_step
        # self.PML_width = int(self.PML_width/self.grid_step) * self.grid_step
        # print(self.PML_width)
        # print(self.grid_step)
        # self.domain_z = int(self.domain_z/self.grid_step) * self.grid_step

        self.cell_size = mp.Vector3(self.domain_x + 2*self.PML_width,
                                        self.domain_y + 2*self.PML_width,
                                        self.domain_z + 2*self.PML_width)

        self.geometry_center = mp.Vector3(0, 0, -(self.cell_size.z/2 - self.z_top_air_gap - self.PML_width - design_specs['d_layers'][-2] - design_specs['d_layers'][-3]/2))

        self.boundary_layers = [mp.PML(self.PML_width)]
        # print( [self.cell_size.x / self.

    def init_sources_and_monitors(self, f, df, use_monitor=True) :

        if df == 0 :
            source = mp.Source(mp.ContinuousSource(f,width=0.1 ),
                               component=mp.Ez,
                               center=mp.Vector3() )
        else :
            source = mp.Source(mp.GaussianSource(f,df),
                               component=mp.Ez,
                               center=mp.Vector3() )

        self.sources.append(source)

        monitor_distance  = self.z_top_air_gap - 0.03

        # nfreq = 1000

        self.monitors = []

        nearfield = mp.Near2FarRegion(mp.Vector3(0, 0, monitor_distance),
                                      size = mp.Vector3(self.domain_x-2*self.grid_step, self.domain_y-2*self.grid_step, 0),
                                      direction = mp.Z)

        # # fluxr = mp.FluxRegion(center=mp.Vector3(0, 0, monitor_distance),
        # #                       size=mp.Vector3(0,0,0),
        # #                       direction=mp.Z)

        self.monitors.append(self.add_near2far(f, 0, 1, nearfield) if use_monitor else None)#, yee_grid=True))

    # def create_blender_primitive(self, scale_factor=1):
    #     blend_name = f"{self.name}.blendtxt"
    #     with open(blend_name,"w") as file:
    #         file.write("//Simulation 00\n")

    #     if obj.__class__ == mp.Block:
    #         name = "cube"

    #         obj_centre = np.array([obj.center.x*scale_factor,
    #                                obj.center.y*scale_factor,
    #                                obj.center.z*scale_factor])
    #         obj_vertices = [np.array([ -.5 * obj.size.x, -.5 * obj.size.y]),
    #                         np.array([ -.5 * obj.size.x, +.5 * obj.size.y]),
    #                         np.array([ +.5 * obj.size.x, +.5 * obj.size.y]),
    #                         np.array([ +.5 * obj.size.x, -.5 * obj.size.y])]


    #         obj_height = obj.size.z

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

            sim_domain = ops.Cube([(self.cell_size.x - self.PML_width) * scale_factor,
                                    (self.cell_size.y - self.PML_width) * scale_factor,
                                    (self.cell_size.z - self.PML_width) * scale_factor], center = True)
            sim_domain = sim_domain.color([.5, .5, .5, .5])
            sim_domain = sim_domain.translate([self.geometry_center.x*scale_factor,
                                                self.geometry_center.y*scale_factor,
                                                self.geometry_center.z*scale_factor])

            # sim_domain.write(scad_name, mode='a')

            # for obj in self.monitors:
            #     center, size = mp.get_center_and_size(obj.where)
            #     plane = ops.Cube([size.x*scale_factor,
            #                       size.y*scale_factor,1e-3*scale_factor], center = True)
            #     plane = plane.color([.5, .5, 0, .5])
            #     plane = plane.translate([center.x*scale_factor,
            #                               center.y*scale_factor,
            #                               center.z*scale_factor])
            #     plane.write(scad_name, mode='a')



#%% geometry and simulation parameters
c0 = 1
wavelength = 0.570
wwidth = .20
f = c0 / wavelength

fmax = c0 / (wavelength - wwidth/2)
fmin = c0 / (wavelength + wwidth/2)
df = 0# fmax - fmin

sim_end = 1e-4

n_eff_l = 1.6642
n_eff_h = 1.7899
n_eff = n_eff_h*.7 + n_eff_l*.3

pattern_type = 'positive'

outcoupler_period = 2*round(wavelength/(n_eff_l+n_eff_h),3)
N_periods = 1
D = 1
charge = 1

t0 = time.time()

# file = 'design_TE_StellaEtAll_2019'
file = 'design_TM_gd3_buriedDBR_onSiO2'

if len(sys.argv) > 1:
    sim_prefix = f"{sys.argv[1]}_"
else:
    sim_prefix = ""

sim_name = f"spiral_outcoupler_{sim_prefix}{file}_{pattern_type}_N{N_periods}_charge{charge}_D{D*1e3:.0f}nm_simend0{sim_end:.1e}"
sim = Simulation(sim_name)

sim.init_geometric_objects(
                multilayer_file = f"Lumerical-Objects/multilayer_design/designs/{file}",
                D = D,
                grating_period = outcoupler_period,
                N_rings = N_periods,
                N_arms  = charge,
                pattern_type=pattern_type)

sim.init_sources_and_monitors(f, df)
mp.verbosity(1)
sim.init_sim()

# sim.create_openscad(scale_factor = 1e3)
# raise ValueError()

date = time.strftime('%y%m%d-%H%M%S')#'211001-121139'#
sim_suffix = f'res{sim.resolution}_{date}'

print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to initiate\n')
#%%
simsize = sim.cell_size
center  = sim.geometry_center


# # # plt.figure(dpi=200)
# fig = plt.figure(dpi=300)
# ax1 = fig.add_subplot(1, 2, 1)
# sim.plot2D( output_plane=mp.Volume(center=center,size=mp.Vector3(0,simsize.y,simsize.z)),
#                 labels=True,
#                 eps_parameters={"interpolation":'none',"cmap":'gnuplot'} )
# ax2 = fig.add_subplot(1, 2, 2)
# sim.plot2D( ax=ax2, output_plane=mp.Volume(size=mp.Vector3(simsize.x,simsize.y)),
#                 labels=True,
#                 eps_parameters={"interpolation":'none',"cmap":'gnuplot'})
# fig.savefig(f'{sim_name}-{sim_suffix}_section.jpg')
# # plt.close()
# plt.show()
# plt.show(block=False)
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

sim.run(mp.at_every(1,print_time),until=10)
# sim.run(until_after_sources=mp.stop_when_fields_decayed(1, mp.Ez, mp.Vector3(), sim_end))

t = np.round(sim.round_time(), 2)

sim.save_near2far(near2far=sim.monitors[0], fname=f'{sim_suffix}_nearfield_t{t}')

ex_near, ey_near = [sim.get_dft_array(sim.monitors[0], field, 0) for field in [mp.Ex, mp.Ey]]

mpo.savemat(f'{sim_name}-{sim_suffix}_nearfield_t{t}.mat', {'Ex': ex_near, 'Ey': ey_near,
                                                            'Lx': sim.monitors[0].regions[0].size.x,
                                                            'Ly': sim.monitors[0].regions[0].size.y})

print(f'\n\nSimulation took {convert_seconds(time.time()-t0)} to run\n')

t1 = time.time()

r = 1e6 #1m
n_freq = 3
res = n_freq/(r/3)

fields = sim.get_farfields(near2far=sim.monitors[0], resolution=res, center=mp.Vector3(0,0,r), size=mp.Vector3(r/3,r/3,0))
Ex, Ey, Ez = [ fields[k] for k in ['Ex', 'Ey', 'Ez']]
mpo.savemat(f'{sim_name}-{sim_suffix}_farfield_t{t}.mat', fields)

print(f'\n\nFar field took {convert_seconds(time.time()-t1)} to compute\n')

del ex_near, ey_near, fields, Ex, Ey, Ez
plt.close()
sim.reset_meep()

sim = Simulation(sim_name)

sim.init_geometric_objects(
                multilayer_file = file,
                D = D,
                grating_period = outcoupler_period,
                N_rings = N_periods,
                N_arms  = charge,
                pattern_type=pattern_type)

sim.init_sources_and_monitors(f, df, use_monitor=False)
mp.verbosity(1)
sim.init_sim()

f = plt.figure(dpi=150)
Animate = mp.Animate2D(sim, fields=mp.Hy, f=f, realtime=False, normalize=True,
                        output_plane=mp.Volume(center=center,size=mp.Vector3(0,simsize.y,simsize.z)))

sim.run(mp.at_every(0.1,Animate),until=10)

filename = f'{sim_name}-{sim_suffix}_section.mp4'
Animate.to_mp4(10,filename)
