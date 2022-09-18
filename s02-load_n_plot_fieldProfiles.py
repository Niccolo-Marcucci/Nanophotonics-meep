#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 18:58:43 2022

@author: ashitaka
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
import os
import sys
import json
import time

files = os.listdir("data")

hashtag = '9772782b76'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file
print()
print(folder)
print()
print()
files = os.listdir(folder)
date = time.strftime('%y%m%d-%H%M%S')

for file in files :
    if file[-4:] == ".mat":
        if file.find("spectra") >= 0:
            data = mpo.loadmat(folder+'/'+file)

#%%
            x = data['field_profile_x'][0]
            y = data['field_profile_y'][0]
            eps = (data['field_profile_Eps'])#)np.flip(np.transpose
            X,Y = np.meshgrid(x,y)
            freqs = 1/ data['field_profile_frequencies'][0]*1e3
            for i in range (len(freqs)):
                Ey = data[f'field_profile_Ey_{i}']
                Ex = data[f'field_profile_Ex_{i}']
                Ez = data[f'field_profile_Ez_{i}']
                fig = plt.figure()
                plt.pcolormesh(X,Y,np.flip(np.transpose(np.log(np.abs(Ez)**2 ))))#+ np.abs(Ex)**2)),axis=0))
                plt.axis('image')
                plt.title(f"{freqs[i]} nm")
                fig.savefig(folder + f'/{file[:-4]}_{date}_mode{i}.png')

            id_x = np.argmin( abs(x))
            id_y = np.argmin( abs(y))
            fig = plt.figure()
            plt.plot(y, np.sqrt(eps[id_x, :]), x, np.sqrt(eps[:, id_y]))
            plt.legend(['x','y'])
            fig = plt.figure()
            plt.pcolormesh(X,Y,np.flip(np.transpose(eps)),cmap='gnuplot')
            plt.axis('image')
            # ax = fig.add_subplot(111, projection='3d')
            # ax.plot_surface(X,Y,(eps),cmap='gnuplot')
