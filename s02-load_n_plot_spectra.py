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

folder = "data/2D_parametric_spacer_220308-174312_301b3cc982"
sim_prefix = "2D_eff_index_cavity_220308-174312_301b3cc982"
period = 280

files = os.listdir(folder)
Ds_list = [[]]
source_positions = []
frst_param = {'name': '_D',
              'list': [[]],
              'label': 'Ds [nm]'}
scnd_param = {'name': '_src',
              'list': [],
              'label': 'Source Position [nm]'}
sim_filelist = [[]]

for file in files :
    if file[:len(sim_prefix)] == sim_prefix and file[-4:] == ".mat":
        if file.find("empty") >= 0:
            data = mpo.loadmat(folder+'/'+file)
            spectrum_empty = data['spectra'][0]
        else :
            suffix = file[len(sim_prefix):]
            first_parameter  = suffix[suffix.find( frst_param['name'] )+2:]
            first_parameter  = int(first_parameter[:first_parameter.find("_")])
            second_parameter = suffix[suffix.find( scnd_param['name'] )+4:]
            second_parameter = int(second_parameter[:second_parameter.find('_')])
            if not any([second_parameter==val for val in scnd_param['list']]):
                scnd_param['list'].append(second_parameter)

                sim_filelist[-1].append(file)
                frst_param['list'][-1].append(first_parameter)
            else:
                j = scnd_param['list'].index(second_parameter)
                sim_filelist[j].append(file)
                frst_param['list'][j].append(first_parameter)

for j, second_parameter in enumerate(scnd_param['list']):

    first_parameter = np.array(frst_param['list'][j])

    sort_idx = first_parameter.argsort(0)
    first_parameter = first_parameter[sort_idx]
    sim_filelist[j] = [sim_filelist[j][idx] for idx in sort_idx]

    image = []

    for k, file in enumerate(sim_filelist[j]) :
        data = mpo.loadmat(folder + '/' + file)
        spectrum = (data['spectra'][0] )
        wavelength = data["wavelength"]
        wavelength = wavelength.reshape((wavelength.size,))
        # sort_idx = wavelength.argsort(0)
        # wavelength = wavelength[sort_idx]
        # spectrumm = spectrum[sort_idx]

        # l=np.int((k-1)/3)
        # # resonance_table = [ [ np.round(1/var[2][l]*1e3, 1), np.int(var[3][l]) ]  for l in range(var[2].size) ]
        # fig = plt.figure(dpi=150,figsize=(10,5))
        # ax = fig.add_subplot(111)
        # legend_str = []
        # ax.plot(wavelength, spectrum)
        # plt.xlim(550,650)
        # # plt.ylim(-2,2)
        # ax.grid(True)
        # plt.xlabel('Wavelength [nm]')
        # plt.ylabel('Transmission')
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
        image.append(np.log10(spectrum))

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

    image = np.array(image)#.transpose()
    fig = plt.pcolor(wavelength, first_parameter , image)
    fig = fig.get_figure()
    fig.set_figheight(6)
    fig.set_figwidth(12)
    fig.set_dpi(150)
    plt.xlabel('wavelength [nm]')
    plt.ylabel(f'{frst_param["label"]}')
    plt.title(f'Period DBR {period}nm, {scnd_param["label"]} {second_parameter}, eff index 1.1')
    fig.savefig(f'{folder}/{sim_prefix}{scnd_param["name"]}{second_parameter}_spacer_dependence_DBRperiod{period}.png')