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
hashtag ='0f198e254a'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file

files = os.listdir(folder)

source_positions = []
scnd_param = {'name': 'point', # '_D',
              'list': [],
              'label': 'effective index'} #'Ds [nm]'}
frst_param = {'name': '_wv', # '_src',
              'list': [],
              'label': 'wavelength (nm)'} # 'Source Position [nm]'}

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
            if file.find("empty") >= 0:
                spectrum_empty = np.mean(abs(data['spectra']),0)
            else :
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
#%%
lambd = np.linspace(.570,.610,500)
inc_sum = np.zeros(lambd.size)

padding = int(5e4);
source_len = int(3e3)
FF_list = []
for j, second_parameter in (enumerate(scnd_param['list'])):
    first_parameter = np.array(frst_param['list'][j])

    # sort first parametr
    sort_idx = first_parameter.argsort(0)
    first_parameter = first_parameter[sort_idx]
    sim_filelist[j] = [sim_filelist[j][idx] for idx in sort_idx]
    if len(first_parameter) > 1:
        print()
    image = []

    # raise
    flag = True
    for k, file in tqdm(enumerate(sim_filelist[j])) :
        data = mpo.loadmat(folder + '/' + file)
        Ez = data['E_z']
        shape = Ez.shape
        if len(shape) < 4 :
            shape = ( *shape, 1,1)
            Ez = Ez.reshape(shape)
        t = np.linspace(0, 400, shape[1]);
        dt = np.diff(t)[0];
        t = np.concatenate( (t, t[-1:] + np.cumsum(dt*np.ones( padding)) ) )
        if flag:
            f = np.linspace(-1/dt/2, 1/dt/2, t.size)
            wv =1/f[(1/f>.55)*(1/f <.63)]
            images = np.zeros( ( len(sim_filelist[j]), wv.size, shape[0]) )
            WV = np.zeros( ( len(sim_filelist[j]), wv.size) )
            WVV = np.zeros( ( len(sim_filelist[j]), wv.size) )
        flag = False
        for i in range(shape[0]-1):#[shape[0]-1]: #  ): # #  [1,5,9,13,16]: # each monitor
            FF = np.zeros(t.size)
            for kk in range(shape[3]):
                 for jj in (range(shape[2])) : # each monitor point
                     # ex = np.concatenate( (E_x[l][i,:,j,k], np.zeros(padding)) )
                     # ex[1:source_len] = 0
                     # ey = np.concatenate( (E_y[l][i,:,j,k], np.zeros(padding)) )
                     # ey[1:source_len] = 0
                     ez = np.concatenate( (Ez[i,:,jj,kk], np.zeros(padding)) )
                     # ez_log = np.log(abs(ez))
                     # ez_log[ez_log==-np.Inf] = 0
                     # p = np.polyfit(t[source_len:shape[1]], ez_log[source_len:shape[1]],1)
                     # ez = ez / np.exp(np.polyval(p, t))
                     ez[:source_len] = 0


                     # wv = np.linspace(.500, .600, 100)
                     # f = 1/wv
                     f = np.linspace(-1/dt/2, 1/dt/2, FF.size)
                     f -= np.diff(f)[0]/2 if np.mod(FF.size,2) else 0 # shift only if the vector is odd
                     FF += abs(np.fft.fftshift(np.fft.fft(ez)))**2
            images[k,:,i] = FF[(1/f>.550)*(1/f <.63)]
        wavelength = 1/f[(1/f>.55)*(1/f <.630)]
        WV[k,:]  = wavelength
        WVV[k,:] = first_parameter[k]*1e-3

    fig = plt.figure()
    image = sum([images[:,:, i] for i in range(shape[0])])
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(WV, WVV, image)
    plt.pcolormesh(WV, WVV, image)
    plt.axis("image")
    # plt.pcolor(wavelength, first_parameter , image)
    plt.plot(WVV[:,0],WVV[:,0])

    fig.set_figheight(6)
    fig.set_figwidth(15)
    fig.set_dpi(100)
    plt.xlabel('wavelength (nm)')
    plt.ylabel(f'{frst_param["label"]}')
    plt.title(f'Period DBR: {period}nm, {scnd_param["label"]}: {second_parameter:.4f}, spacer: 640nm')
    # fig.savefig(folder + f'/maps/sim_2D_{date}_{scnd_param["name"]}{second_parameter+360:.0f}_intensity_map.png')

    # plt.close()

    spectra =  np.zeros( (shape[0]-1 , lambd.size) )

    fig = plt.figure()
    fig.set_figheight(6)
    fig.set_figwidth(12)
    def run_parallel(image):
        return itp.griddata((WV.reshape(WV.size), WVV.reshape(WV.size)), image.reshape(WV.size), (lambd, lambd))
    with Pool(1) as parfor:
        output = parfor.map(run_parallel, (images[:,:,i] for i in range(shape[0])))#range(len(spectra[:,0])) ))
    # output = [run_parallel(images[:,:,i]) for i in ]
    # for i in  tqdm(range(0,len(data['spectra']))):
    #     spectra[i,:] = itp.griddata((WV.reshape(WV.size), WVV.reshape(WV.size)), images[:,:,i].reshape(WV.size), (lambd, lambd))
    spectra = np.array(output)
    # io.savemat(folder + f'/mat_files/sim_2D_{date}_{scnd_param["name"]}{second_parameter+360:.0f}_opposite_monitor_spectrum.mat',
    #             {"diagonal_wavelength":lambd, "diagonal_intensity": spectra[0], "wv_n_eff": WV, "wavelength" : WVV, "map": image})
    intensity =  sum([spectra[i,:] for i in range(shape[0]-1,)])#[3,7,11,15]])
    plt.plot(lambd, intensity )
    plt.xlabel('wavelength (nm)')
    plt.ylabel('intensity (a.u.)')
    plt.grid(True)
    # images.append(image)
    plt.title(f'Period DBR: {period}nm, source_{scnd_param["label"]}: {second_parameter:.0f}, spacer: 560nm')
    # fig.savefig(folder + f'/intensity/sim_2D_{date}_{scnd_param["name"]}{second_parameter+360:.0f}_intensity.png')
    plt.close()
    inc_sum += intensity
#%%
fig = plt.figure()
fig.set_figheight(6)
fig.set_figwidth(12)
plt.plot(lambd, inc_sum )
plt.xlabel('wavelength (nm)')
plt.ylabel('intensity (a.u.)')
plt.grid(True)
# images.append(image)
plt.title(f'Period DBR: {period}nm, spacer: 560nm')
# fig.savefig(folder + f'/sim_2D_{date}_incoerent_sum_intensity.png')
# fig = plt.figure()
# image = np.array(images[0] + images[1])#.transpose()
# plt.pcolor(wavelength, first_parameter , image)

# fig.set_figheight(3)
# fig.set_figwidth(12)
# fig.set_dpi(100)
# plt.xlabel('wavelength [nm]')
# plt.ylabel(f'{frst_param["label"]}')
# plt.title(f'Period DBR {period}nm, {scnd_param["label"]} {second_parameter}, eff index 1.1')