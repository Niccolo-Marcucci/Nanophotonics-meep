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
from tqdm import tqdm

def ft(x, t, f):
    dt =  (t[1]-t[0])

    # fourier transform definition
    # F = np.zeros(f.size, dtype=complex)
    # for i in range(f.size) :
    #     F[i] = sum( x * np.exp( -1j*2*np.pi * f[i] * t) ) *dt

    # expoit matrix product to speed up
    f = np.matrix(f)
    t = np.matrix(t)
    x = np.matrix(x)
    F = np.array(np.exp( -1j*2*np.pi * f.transpose() * t) * x.transpose()) * dt
    return F.reshape((F.size,))

files = os.listdir("data")

hashtag ='42dd1ce876'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file
print()
print(folder)
print()
print()
files = os.listdir(folder)
E_x = []
E_y = []
E_z = []
titles = []
i = 0
for file in files :
    if file[-4:] == ".mat":
        if file.find("spectra") >= 0:
            i += 1
            print(file)
            data = mpo.loadmat(folder+'/'+file)
            if file.find("empty") >= 0:
                spectrum_empty = data['spectra']
                FT_x_empty = data['E_x']
                FT_y_empty = data['E_y']
            else :
                titles.append(file)
                E_x.append(data['E_x'])
                E_y.append(data['E_y'])
                E_z.append(data['E_z'])


if i == 0 :
    raise FileNotFoundError(f"No spectra file was found for hash {hashtag}")

#%%
fig = plt.figure(dpi=150,figsize=(10,6))
ax1 = fig.add_subplot(111)
plt.title("FFT ")#, monitor a 0-90")#file
fig2 = plt.figure(dpi=150,figsize=(10,6))
ax2 = fig2.add_subplot(111)

#%%
flag = True
padding = int(0);
source_len = int(12e3)
FF_list = []
for l in range(len( E_x[:])):

    shape = E_z[l].shape
    if len(shape) < 4 :
        shape = ( *shape, 1,1)
        E_x[l] = E_x[l].reshape(shape)
        E_y[l] = E_y[l].reshape(shape)
        E_z[l] = E_z[l].reshape(shape)

    if flag :
        t = np.linspace(0, 400, shape[1]);
        dt = np.diff(t)[0];
        t = np.concatenate( (t, t[-1:] + np.cumsum(dt*np.ones( padding)) ) )
        FF = np.zeros(t.size)
        field = np.zeros( (1000, shape[2],shape[3]) )
    flag = False

    for i in [shape[0]-1]: # tqdm(range(shape[0]-1)): # #  [1,5,9,13,16]: # each monitor
        for j in tqdm(range(shape[2])) : # each monitor point
            for k in range(shape[3]):
                # ex = np.concatenate( (E_x[l][i,:,j,k], np.zeros(padding)) )
                # ex[1:source_len] = 0
                # ey = np.concatenate( (E_y[l][i,:,j,k], np.zeros(padding)) )
                # ey[1:source_len] = 0
                ez = np.concatenate( (E_z[l][i,:,j,k], np.zeros(padding)) )
                ez_log = np.log(abs(ez))
                ez_log[ez_log==-np.Inf] = 0
                p = np.polyfit(t[source_len:shape[1]], ez_log[source_len:shape[1]],1)
                ez = ez / np.exp(np.polyval(p, t))
                ez[:source_len] = 0


                # wv = np.linspace(.500, .600, 100)
                # f = 1/wv
                f = np.linspace(-1/dt/2, 1/dt/2, FF.size)
                f -= np.diff(f)[0]/2 if np.mod(FF.size,2) else 0 # shift only if the vector is odd
                FF = abs(np.fft.fftshift(np.fft.fft(ez)))**2 #+ abs(np.fft.fftshift(np.fft.fft(ex)))**2#
                # FF = abs(ft(ez, t, f))**2 #+ abs(np.fft.fftshift(np.fft.fft(ex)))**2#
                FF_list.append(FF[(1/f>.57)*(1/f <.62)])

                wvs = 1/f[(1/f>.57)*(1/f <.62)]

                for s, wv in enumerate(wvs):
                    idx  = np.argmin(abs(1/f - wv))
                    field[s,j,k] = FF[idx]
                # ax2.semilogy(t, abs(ez)**2) #+abs(ey)**2)


        ax1.plot(1/f*1e3, FF, '-' if i<10 else '--')

plt.draw()
#%%
for ax in [ax1, ax2]:
    ax.grid(True)
    ax.minorticks_on()
    # ax.grid(True, 'minor','x')
ax1.set_xlim(550, 620)
ax1.set_xlabel('Wavelength [nm]')
ax2.set_xlabel('time (s)')
ax1.set_ylabel('Transmission')
ax2.set_ylabel('field')
ax1.legend(range(shape[0]-1))
date = time.strftime('%y%m%d-%H%M%S')

fig.savefig(folder + '/' + date + '_fft.png')

raise
# plt.legend()
#%%
# for s, wv in enumerate(wvs): #range(100,200,10) :
#     # plt.figure()
#     # plt.imshow(field[s,:,:])
#     plt.imsave(folder+f"/fields/field_{wv*1e3:.1f}.png", field[s,:,:])
#%%
# FF_list.append(FF_list[0])
im = np.transpose(np.array(FF_list))# [(1/f>.57)*(1/f <.595), :]

wv = 1/ f[(1/f>.57)*(1/f <.62 )] - .57

theta = np.linspace(-np.pi/2,np.pi/2,9)[1:] # - np.pi/32
THT, WV = np.meshgrid(theta, wv)
x = WV * np.cos(THT)
y = WV * np.sin(THT)

plt.figure()
plt.pcolormesh(x, y, im)
plt.axis('image')