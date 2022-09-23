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
hashtag ='9388d33171'#'e20d2ea866'#5cab4c8862'#'7bae1ab6b6'#''62ef19ee4d'#'4dc3971d95'#'8a593f9138'#'2cb6bfb1fa'

for file in files :
    if file.find( hashtag ) >= 0:
        folder = 'data/' + file

files = os.listdir(folder)

source_positions = []
scnd_param = {'name': '_angle', # '_D',
              'list': [],
              'label': 'angle (°)'} #'Ds [nm]'}
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
second_parameter = np.array([-angle/np.pi**2 if angle <= 0 else 180-angle/np.pi**2 for angle in scnd_param['list']]) # valid only for angle!!!

if second_parameter.size > 1:
    sort_idx = second_parameter.argsort(0)
    sim_filelist = [sim_filelist[idx] for idx in sort_idx]
    scnd_param['list'] = [round(second_parameter[idx]) for idx in sort_idx]
    frst_param['list'] = [frst_param['list'][idx] for idx in sort_idx]
# fig = plt.figure(dpi=150,figsize=(10,5))
# ax = fig.add_subplot(111)
#%%

lambd = np.linspace(570,610,500)
inc_sum = np.zeros(lambd.size)
flag = True
for j, second_parameter in tqdm(enumerate(scnd_param['list'])):
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
    images = np.zeros( ( len(sim_filelist[j]), spectrum_empty.size, len(data['spectra'])) )
    WV = np.zeros( ( len(sim_filelist[j]), spectrum_empty.size) )
    WVV = np.zeros( ( len(sim_filelist[j]), spectrum_empty.size) )
    for k, file in enumerate(sim_filelist[j]) :
        data = mpo.loadmat(folder + '/' + file)
        # fig = plt.figure(dpi=150,figsize=(10,5))
        # ax = fig.add_subplot(111)
        for i in [0] : #range(len(data['spectra'])):# [0]:#[7,11,15,] : # [7,15] : ##
            images[k,:,i] = (abs(data['spectra'][i])**2) # np.abs(data['FT_x'])**2 + np.abs(data['FT_y'])**2 #
        wavelength = data["wavelength"][0]
        WV[k,:]  = wavelength
        WVV[k,:] = 1 / ( (1/wavelength[-1] + 1/wavelength[1])/2 )
        # ax.plot(wavelength, np.log(np.abs(spectrum)))
        # plt.legend(range(0,16))

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
    image = sum([images[:,:, i] for i in range(len(data['spectra']))])
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
    plt.title(f'Period DBR: {period}nm, source_{scnd_param["label"]}: {second_parameter:.0f}, spacer: 560nm')
    fig.savefig(folder + f'/sim1D_{date}_{scnd_param["name"]}{second_parameter:.0f}_intensity_map.png')

    plt.close()
    if flag :
        image_sum = image
        flag = False
    else:
        image_sum += image
    spectra =  np.zeros( ( len(data['spectra']), lambd.size) )

    fig = plt.figure()
    fig.set_figheight(6)
    fig.set_figwidth(12)
    def run_parallel(image):
        return itp.griddata((WV.reshape(WV.size), WVV.reshape(WV.size)), image.reshape(WV.size), (lambd, lambd))
    with Pool(6) as parfor:
        output = parfor.map(run_parallel, (images[:,:,i] for i in range(len(data['spectra'])) )) # [0]))#
    # output = [run_parallel(images[:,:,0])]
    # for i in  tqdm(range(0,len(data['spectra']))):
    #     spectra[i,:] = itp.griddata((WV.reshape(WV.size), WVV.reshape(WV.size)), images[:,:,i].reshape(WV.size), (lambd, lambd))
    spectra = np.array(output)
    # io.savemat(folder + f'/mat_files/sim_2D_{date}_{scnd_param["name"]}{second_parameter:.0f}_opposite_monitor_spectrum.mat',
    #             {"diagonal_wavelength":lambd, "diagonal_intensity": spectra, "monitors_angles": np.linspace(-180,180,17)[1:], "wv_n_eff": WV, "wavelength" : WVV, "maps": images})
    intensity =  sum([spectra[i,:] for i in range(len(data['spectra']))])#[3,7,11,15]])
    plt.plot(lambd, intensity )
    plt.xlabel('wavelength (nm)')
    plt.ylabel('intensity (a.u.)')
    plt.grid(True)
    # images.append(image)
    plt.title(f'Period DBR: {period}nm, source_{scnd_param["label"]}: {second_parameter:.0f}, spacer: 560nm')
    fig.savefig(folder + f'sim1D_{date}_{scnd_param["name"]}{second_parameter:.0f}_intensity.png')
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

fig.savefig(folder + f'/ortog_monitors_sim_2D_{date}_incoerent_sum_intensity.png')

fig = plt.figure()
plt.pcolormesh(WV, WVV, image_sum)

fig.set_figheight(6)
fig.set_figwidth(15)
fig.set_dpi(100)
plt.xlabel('wavelength (nm)')
plt.ylabel(f'{frst_param["label"]}')
plt.title(f'Period DBR: {period}nm, spacer: 560nm')
fig.savefig(folder + f'/ortog_monitors_sim_2D_{date}_incoerent_sum_intensity_map.png')