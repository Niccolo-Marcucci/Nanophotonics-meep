#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 11:21:57 2022

@author:
"""

import meep as mp
import numpy as np
# import h5py as h5
#import scipy as sp
from scipy import optimize as opt
from scipy import interpolate as itp, io
from matplotlib import pyplot as plt
from multiprocessing import Pool
# from mpl_toolkits.mplot3d import Axes3D
import meep_objects as mpo
import os
import sys
import json
import time
from tqdm import tqdm
date = time.strftime('%y%m%d-%H%M%S')


files = os.listdir("data")
hashtag ='68fa746847'#'1c62b710f2'#'038dd156ce'#'7edc4aea0b'#4fb109da14'#'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file

files = os.listdir(folder)

source_positions = []
scnd_param = {'name': '_n_eff_h', # '_D',
              'list': [],
              'label': 'angle (°)'} #'Ds [nm]'}
frst_param = {'name': '_n_eff_l', # '_src',
              'list': [],
              'label': 'n_eff wavelength (nm)'} # 'Source Position [nm]'}

# sim_filelist will be a two dimensional list containing all useful files
# (same as frst_param["list"]). In other words, each element of the list is a
# list containing all files with referred to a given value of scnd_param
sim_filelist = []
spectrum_empty = []
i = 0
for file in files :
    if file[-4:] == ".mat":
        if file.find("spectra") >= 0:
            i += 1
            data = mpo.loadmat(folder+'/'+file)

            suffix = file[file.find(hashtag)+len(hashtag):]
            first_parameter  = suffix[suffix.find( frst_param['name'] ) + len(frst_param['name']):]
            first_parameter  = float(first_parameter[:first_parameter.find("_")])
            second_parameter = suffix[suffix.find( scnd_param['name'] ) + len(scnd_param['name']):]
            second_parameter = float(second_parameter[:second_parameter.find('_')])
            if not any([second_parameter==val for val in scnd_param['list']]):
                # if the same value of second_parameter does not exist creat a new entry
                scnd_param['list'].append(second_parameter)

                sim_filelist.append([file])
                frst_param['list'].append([first_parameter])
            else:
                j = scnd_param['list'].index(second_parameter)
                sim_filelist[j].append(file)
                frst_param['list'][j].append(first_parameter)

if i == 0 :
    raise FileNotFoundError(f"No spectra file was found for hash {hashtag}")


period = 280

#%% sort second parameter list
second_parameter = np.array(scnd_param['list'])
if second_parameter.size > 1:
    sort_idx = second_parameter.argsort(0)
    sim_filelist = [sim_filelist[idx] for idx in sort_idx]
    scnd_param['list'] = [scnd_param['list'][idx] for idx in sort_idx]
    frst_param['list'] = [frst_param['list'][idx] for idx in sort_idx]
# fig = plt.figure(dpi=150,figsize=(10,5))
# ax = fig.add_subplot(111)

resonances =np.zeros( (len(scnd_param['list']), len(scnd_param['list'])) )
qualities = np.zeros( (len(scnd_param['list']), len(scnd_param['list'])) )

for j, second_parameter in tqdm(enumerate(scnd_param['list'])):
    first_parameter = np.array(frst_param['list'][j])

    # sort first parametr
    sort_idx = first_parameter.argsort(0)
    first_parameter = first_parameter[sort_idx]
    sim_filelist[j] = [sim_filelist[j][idx] for idx in sort_idx]

    for k, file in enumerate(sim_filelist[j]) :
        data = mpo.loadmat(folder + '/' + file)
        resonances[j,k] = data["resonance_table_t170.0"][0,0]
        qualities[j,k] = data["resonance_table_t170.0"][0,1]

d = np.linspace(0, 65, 30)
X, Y = np.meshgrid(scnd_param['list'], first_parameter)
X, Y = np.meshgrid(d[:-1], d[1:], indexing='xy')

resonances[resonances==0.] = np.NaN
# resonances[Y>51] = np.NaN
#%%
def get_contour_verts(cn):
    lines = []
    # for each contour line
    for cc in cn.collections:
        # for each separate section of the contour line
        for pp in cc.get_paths():
            xy = []
            # for each segment of that section
            for vv in pp.iter_segments():
                xy.append(vv[0])
            # exclude short segments
            if len(xy) > 10:
                lines.append( np.stack(xy))

    return np.array(lines)

plt.close('all')

levels_wv = [584, 596 ]
levels_q = np.linspace(2,2.8, 11)

fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, resonances)
plt.pcolormesh(X, Y, resonances)
plt.colorbar()
p = plt.contour(X, Y, resonances, levels_wv, cmap='hot')
plt.title("resonance wavelength")
plt.xlabel("lower thickness nm")
plt.ylabel("upper thickness nm")
# plt.x  lim([min(X.min(), Y.min()), max(X.max(), Y.max())])
# plt.ylim([min(X.min(), Y.min()), max(X.max(), Y.max())])
lines_wv = get_contour_verts(p)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
plt.pcolormesh(X, Y, np.log10(qualities))
plt.colorbar()
p = plt.contour(X, Y,  np.log10(qualities), levels_q, cmap='hot')
plt.title("resonance quality")
plt.xlabel("lower thickness nm")
plt.ylabel("upper thickness nm")
lines_q = get_contour_verts(p)

a = [ ax1.plot(line[:,0], line[:,1], 'r' ) for line in lines_q]
b = [ ax2.plot(line[:,0], line[:,1], 'r' ) for line in lines_wv]

# plt.figure(fig), plt.savefig(folder + "/resonances_wavelengths_map.png")
# plt.figure(fig2), plt.savefig(folder + "/resonances_qualities_map.png")
# ax1.plot([2,15],[31, 40], '+k', linewidth=3)

points = []
for line_q in lines_q:
    line_q_ = itp.interp1d(line_q[:,0],line_q[:,1])
    for line_wv in lines_wv:
        x_min = max(min(line_wv[:,0]), min(line_q[:,0]))
        x_max = min(max(line_wv[:,0]), max(line_q[:,0]))

        line_wv_ = itp.interp1d(line_wv[:,0],line_wv[:,1])

        f = lambda x: line_wv_(x) - line_q_(x)

        try:
            root = opt.root_scalar(f, x0=10, x1 = 11, bracket=[x_min, x_max]).root
        except ValueError:
            root = np.NaN
        else:
            points.append([root, line_wv_(root).item()] )
            ax1.plot(points[-1][0], points[-1][1], '+k')
            ax2.plot(points[-1][0], points[-1][1], '+k')
points = np.array(points)
print(points)
# plt.figure()
# for i in range(len(levels_wv)):
#     for k in range(2):
#         plt.plot(np.linspace(2,3, 11),points[:,i,k])



# io.savemat(folder + f"/cross_points_{hashtag}.mat", {"points": points, "levels_wv": levels_wv, "levels_q": levels_q})
# io.savemat(f"cross_points_{hashtag}.mat", {"points": points, "levels_wv": levels_wv, "levels_q": levels_q})
