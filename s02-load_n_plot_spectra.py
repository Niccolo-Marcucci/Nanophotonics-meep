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
date = time.strftime('%y%m%d-%H%M%S')


files = os.listdir("data")
hashtag ='e5dff2b4e3'#'5cab4c8862'#'e20d2ea866'#'62ef19ee4d'#'4dc3971d95'#'8a593f9138'#'2cb6bfb1fa'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file

files = os.listdir(folder)

source_positions = []
scnd_param = {'name': '_D', # '_D',
              'list': [],
              'label': 'specer (nm)'} #'Ds [nm]'}
frst_param = {'name': '_wv', # '_src',
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
            wavelength = data["wavelength"][0]
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

if len(spectrum_empty) == 0:
    spectrum_empty = np.ones((wavelength.size,))

period = 280

#%% sort second parameter list
second_parameter = np.array(scnd_param['list'])
if second_parameter.size > 1:
    sort_idx = second_parameter.argsort(0)
    sim_filelist = [sim_filelist[idx] for idx in sort_idx]
    scnd_param['list'] = [scnd_param['list'][idx] for idx in sort_idx]
# fig = plt.figure(dpi=150,figsize=(10,5))
# ax = fig.add_subplot(111)
#%%
images = []
for j, second_parameter in enumerate(scnd_param['list']):
    first_parameter = np.array(frst_param['list'][j])

    # sort first parametr
    sort_idx = first_parameter.argsort(0)
    first_parameter = first_parameter[sort_idx]
    sim_filelist[j] = [sim_filelist[j][idx] for idx in sort_idx]

    image = []

    # spectrum = np.ones(len(sim_filelist[j]))
    # wavelength = np.ones(len(sim_filelist[j]))
    # for k, file in enumerate(sim_filelist[j]) :
    #     data = mpo.loadmat(folder + '/' + file)

    #     sp = 0
    #     FT_x = 0
    #     FT_y = 0
    #     for i in [16] : #range(len(spectrum)): # [7,15] : ##
    #         sp += abs(data['spectra'][0])
    #         FT_x += max(data['E_x'][0])
    #         FT_y += max(data['E_y'][0])

    #     sp = (np.abs(FT_x)**2 + np.abs(FT_y)**2)
    #     spectrum[k] = ( sp /spectrum_empty)
    #     wavelength[k] = data["wavelength"][0]

    # legend_str = []
    # ax.plot(wavelength, np.log(np.abs(spectrum)), '.')
    # plt.xlim(550,650)
    # # plt.ylim(-2,2)
    # ax.grid(True)
    # plt.xlabel('Wavelength [nm]')
    # plt.ylabel('Transmission')
    # date = time.strftime('%y%m%d-%H%M%S')
    # fig.savefig(folder + '/' + date + '_spectrum.png')

    # raise
    image = np.zeros( ( len(sim_filelist[j]), 1000) )
    WV = np.zeros( ( len(sim_filelist[j]), 1000) )
    WVV = np.zeros( ( len(sim_filelist[j]), 1000) )
    for k, file in enumerate(sim_filelist[j]) :
        data = mpo.loadmat(folder + '/' + file)
        spectrum = abs(data['spectra'][3] )/spectrum_empty #can devide by spectrum_emty even i have only one
        # fig = plt.figure(dpi=150,figsize=(10,5))
        # ax = fig.add_subplot(111)
        for i in [3, 7,11,15,] : #  range(0,len(data['spectra'])):#[7,15] : ##
            spectrum += abs(data['spectra'][i])/spectrum_empty #can devide by spectrum_emty even i have only one
        wavelength = data["wavelength"][0]
        # ax.plot(wavelength, np.log(np.abs(spectrum)))
        # plt.legend(range(0,16))
        # wavelength = wavelength.reshape((wavelength.size,))
        # sort_idx = wavelength.argsort(0)
        # wavelength = wavelength[sort_idx]
        # spectrumm = spectrum[sort_idx]
        WV[k,:]  = wavelength
        WVV[k,:] = 1 / ( (1/wavelength[-1] + 1/wavelength[1])/2 )
        image[k,:] = spectrum

        # l=np.int((k-1)/3)
        # # resonance_table = [ [ np.round(1/var[2][l]*1e3, 1), np.int(var[3][l]) ]  for l in range(var[2].size) ]
        #%%
        # fig = plt.figure(dpi=150,figsize=(10,5))
        # ax = fig.add_subplot(111)
        # legend_str = []
        # ax.plot(wavelength, np.log(np.abs(spectrum)), '-')
        # plt.xlim(550,650)
        # # plt.ylim(-2,2)
        # ax.grid(True)
        # plt.xlabel('Wavelength [nm]')
        # plt.ylabel('Transmission')
        #%%
        # plt.title(f'D = {D}')
        # legend_str.append(f"sourcePos {tuple_list[k][6]*1e3:.0f}")
        # plt.title(f'n_eff_h={tuple_list[k][1]:.2f};   DBR_period={tuple_list[k][4]*1e3:.0f};   D={tuple_list[k][3]/tuple_list[k][2]:.2f}*DBR_period')
        # # ax2 = fig.add_subplot(336)
        # # # plt.title('Table of the resonances')
        # # collabel=[ "Wavelength", "Quality"]
        # # rowlabel=[ f'{i}' for i,_ in enumerate(var["resonance_table_t300.0"])]#[ f'({np.int(np.rad2deg(angles[i][0]))},{np.int(np.rad2deg(angles[i][1]))})' for i in indeces]
        # # ax2.axis('tight')
        # # ax2.axis('off')
        # # the_table = ax2.table(cellText=var["resonance_table_t300.0"], colLabels=collabel, rowLabels=rowlabel,loc='center')
        # plt.legend(legend_str)

        # # plt.show()
        # fig.savefig(f'{names[k]}_spectrum.png')
        # plt.pause(1)
        # plt.close(fig)

        # fig = plt.figure(10*k,dpi=150,figsize=(10,5))
        # ax = fig.add_subplot(111)
        # ax.plot(var["field_profile_x"], var["field_profile_Eps"])
        # field = var["field_profile_Ey_0"]
        # field = np.abs(field/np.max(field))**2
        # ax.plot(var["field_profile_x"], field)
        # ax.grid(True)
        # plt.xlabel('x [um]')
        # plt.legend(["dielectric_constant", "field_profile"])
        # plt.title(f'n_eff_h={tuple_list[k][1]:.2f};   DBR_period={tuple_list[k][4]*1e3:.0f}; D={tuple_list[k][3]/tuple_list[k][4]:.2f}*DBR_period')
        # fig.savefig(f'{names[k]}_spectrumfield_profile.png')
    fig = plt.figure()
    # image = np.array(image)#.transpose()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(WV, WVV, image)
    plt.pcolormesh(WV, WVV, image)
    # plt.pcolor(wavelength, first_parameter , image)
    plt.plot(WVV[:,0],WVV[:,0])

    fig.set_figheight(6)
    fig.set_figwidth(12)
    fig.set_dpi(100)
    plt.xlabel('wavelength (nm)')
    plt.ylabel(f'{frst_param["label"]}')
    plt.title(f'Period DBR {period}nm, {scnd_param["label"]} {second_parameter}, eff index 1.1')
    fig.savefig(folder + f'/{date}_{scnd_param["name"]}{second_parameter}_intensity_map.png')

    fig = plt.figure()
    fig.set_figheight(6)
    fig.set_figwidth(12)
    lambd = np.linspace(570,610,1000)
    intensity = itp.griddata((WV.reshape(WV.size), WVV.reshape(WV.size)), image.reshape(WV.size), (lambd, lambd))
    plt.plot(lambd, intensity )
    plt.xlabel('wavelength (nm)')
    plt.ylabel('intensity (a.u.)')
    plt.grid(True)
    # images.append(image)
    fig.savefig(folder + f'/{date}_{scnd_param["name"]}{second_parameter}_intensity.png')
#%%
# fig = plt.figure()
# image = np.array(images[0] + images[1])#.transpose()
# plt.pcolor(wavelength, first_parameter , image)

# fig.set_figheight(3)
# fig.set_figwidth(12)
# fig.set_dpi(100)
# plt.xlabel('wavelength [nm]')
# plt.ylabel(f'{frst_param["label"]}')
# plt.title(f'Period DBR {period}nm, {scnd_param["label"]} {second_parameter}, eff index 1.1')