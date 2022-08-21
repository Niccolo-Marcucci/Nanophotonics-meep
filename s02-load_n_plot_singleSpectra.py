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

hashtag = '17e3490176'#8218750129'#'b807d7cbe0'#'e20d2ea866'#'d3dc776849'#'907626023e' #37d4a928c2'# '3bdc9c5579'#'d8d203aecf' #''a919609ed6'#'4b4005295f' #

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file
print()
print(folder)
print()
print()
files = os.listdir(folder)

spectrum_empty = np.zeros(1)
i = 0
for file in files :
    if file[-4:] == ".mat":
        if file.find("spectra") >= 0:
            i += 1
            data = mpo.loadmat(folder+'/'+file)
            wavelength = data["wavelength"][0]
            if file.find("empty") >= 0:
                spectrum_empty = data['spectra'][0]
            else :
                spectrum = data['spectra'][0]

if i == 0 :
    raise FileNotFoundError(f"No spectra file was found for hash {hashtag}")
elif spectrum_empty.size == 1:
    spectrum_empty = np.ones((wavelength.size,))
period = 280

#%%
fig = plt.figure(dpi=150,figsize=(10,5))
ax = fig.add_subplot(211)
ax.plot(wavelength, spectrum, wavelength, spectrum_empty)
# plt.xlim(550,650)
# plt.ylim(-2,2)
ax.grid(True)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmission (a.u.)')
ax2 = fig.add_subplot(212)
ax2.plot(wavelength, np.log(spectrum/spectrum_empty))
# plt.xlim(550,650)
# plt.ylim(-2,2)
ax2.grid(True)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Transmission')

date = time.strftime('%y%m%d-%H%M%S')
fig.savefig(folder + '/' + date + '_spectrum.png')
#%%