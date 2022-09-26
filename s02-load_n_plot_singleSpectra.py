#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 16:17:30 2022

@author:
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

hashtag ='eab3b1d422'#'ac7da1b8e1'#'970756bcfb'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file
print()
print(folder)
print()
print()
files = os.listdir(folder)
files.sort()

lista_spectra = []
titles = []
spectrum_empty = np.zeros(1)
i = 0
for file in files :
    if file[-4:] == ".mat":
        if file.find("spectra") >= 0:
            i += 1
            print(file)
            data = mpo.loadmat(folder+'/'+file)
            wavelength = data["wavelength"][0]
            if file.find("empty") >= 0:
                spectrum_empty = data['spectra']
                FT_x_empty = data['FT_x'][0]
                FT_y_empty = data['FT_y'][0]
            else :
                titles.append(file)
                spectrum = data['spectra']
                lista_spectra.append(data['spectra'])
                # FT_x = data['FT_x'][0]
                # FT_y = data['FT_y'][0]

if i == 0 :
    raise FileNotFoundError(f"No spectra file was found for hash {hashtag}")

if spectrum_empty.size == 1:
    spectrum_empty = np.ones((wavelength.size,))
period = 280

#%%
fig = plt.figure(dpi=150,figsize=(10,6))
ax1 = fig.add_subplot(111)
plt.title("dipolo lungo XY")#, monitor a 0-90")#file
# ax2 = fig.add_subplot(212)

sp0 = np.zeros(wavelength.shape)
# for i in range(len(spectrum_empty)):#
#     sp0 += abs(spectrum_empty[i])/len(spectrum_empty)
#%%
# sp_centre =  np.abs(FT_x)**2 + np.abs(FT_y)**2
for spectrum in lista_spectra[:]:
    sp = np.zeros(wavelength.shape)
    for i in [0,2] : #range(len(spectrum))[:]: # [[7] : # [7,15] : ##
        sp = (abs(spectrum[i]))#/spectrum_empty[i]/len(spectrum)

        # sp = sp/max(sp)
        # sp_centre = sp_centre/max(sp_centre)
        ax1.plot(wavelength, sp, '-' if i < 10 else '--')

        # ax2.plot(wavelength, np.log(sp), '-' if i < 10 else '--')
        plt.draw()
    # ax1.plot(wavelength, sp_centre)
    # ax2.plot(wavelength, np.log(sp_centre))

#%%
for ax in [ax1, ]:
    ax.grid(True)
    ax.minorticks_on()
    ax.grid(True, 'minor','x')
    ax.set_xlim(540, 640)
    ax.set_xlabel('Wavelength [nm]')
ax1.set_ylabel('Transmission')
# ax2.set_ylabel('Transmission - log scale')
date = time.strftime('%y%m%d-%H%M%S')
plt.axes(ax1), plt.legend([f"source {i}° - monitor {j}°" for i in [0, 90] for j in [0, 90]])
fig.savefig(folder + '/' + date + '_spectrum.png')

# plt.legend()