#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:46:00 2021

@author: ashitaka
"""


import meep as mp
import numpy as np
#import scipy as sp
#from scipy import optimize as op
from scipy import interpolate as itp
from matplotlib import pyplot as plt
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D
import meep_objects as mpo
import io
import sys
import time
import h5py as h5


dataset = h5.File('spiral_outcoupler_design_TM_gd3_buriedDBR_onSiO2_positive_N5_charge1_D5000nm_simend01.0e-04-res66_220111-101700_nearfield_t30.0.h5','r')
# field = [np.array(dataset.get(key)) for key in dataset.keys()]
# ex_near = np.array(dataset.get("ex0.r")) + 1j*np.array(dataset.get("ex0.i"))
# ey_near = np.array(dataset.get("ey0.r")) + 1j*np.array(dataset.get("ey0.i"))

# fields = mpo.loadmat("polSplitter_design_TM_gd3_buriedDBR_onSiO2_N9_Dphi60_sigma-1-res44_211002-124621_nearfield_t9.045454978942871")
# ex_near, ex_near = [ fields[k] for k in ['Ex', 'Ey']]

fields = mpo.loadmat("polSplitter_design_TM_gd3_buriedDBR_onSiO2_positive_N9_Dphi60_sigma-1-res66_220105-132026_nearfield_t9.05_MAT.mat")

# Ex, Ey, Ez = [ fields[k] for k in ['Ex', 'Ey', 'Ez']]
# ex_near, ey_near = [ fields[k] for k in ['Ex', 'Ey']]
# Ex, Ey = [ fields[k] for k in ['Ex', 'Ey']]
Ex, Ey = [ fields[k] for k in ['Ex_far', 'Ey_far']]

# Ex = Ex[:,1:]
# Ey = Ey[1:,:]
r = 1e6 #1m
n_freq = 150
res = n_freq/r

x_far = np.linspace(-r, r, Ex.shape[0])
x_far = x_far / np.sqrt(x_far**2 + r**2)
y_far = np.linspace(-r, r, Ey.shape[0])
y_far = y_far / np.sqrt(y_far**2 + r**2)
x_far = np.linspace(-1, 1, Ex.shape[0])
y_far = np.linspace(-1, 1, Ey.shape[0])

X_far, Y_far = np.meshgrid(x_far, y_far)
r = np.sqrt(X_far**2 + Y_far**2)
theta = np.arcsin(r)
phi = np.arctan2(X_far, Y_far)

ex_far = Ex#**2 + (-Ez * np.cos(phi))**2)
ey_far = Ey#np.sqrt(Ey**2 + (-Ez * np.sin(phi))**2)

er = np.sqrt(2)/2 * ex_far + np.sqrt(2)/2 * ey_far * np.exp(-1j*np.pi/2);
el = np.sqrt(2)/2 * ex_far + np.sqrt(2)/2 * ey_far * np.exp(+1j*np.pi/2);

S3 = 1j*(ex_far * np.conj(ey_far) - ey_far * np.conj(ex_far));
S0 = (np.abs(ex_far)**2 + np.abs(ey_far)**2);
chi = 0.5*np.arcsin( np.real(S3)/S0);

# E = sqrt(real(Ex)^2+real(Ey)^2)+1i*sqrt(imag(Ex).^2+imag(Ey).^2);


mpo.plot_image(x_far, y_far, S0)# + np.abs(Ez)**2))

mpo.plot_image(x_far, y_far, (np.abs(er)**2), cmap='hot',vmax=S0.max(),vmin=0)
mpo.plot_image(x_far, y_far, (np.abs(el)**2), cmap='hot',vmax=S0.max(),vmin=0)
mpo.plot_image(x_far, y_far, (np.angle(er)), cmap='hot')
mpo.plot_image(x_far, y_far, (np.angle(el)), cmap='hot')
mpo.plot_image(x_far, y_far, (np.real(np.tan(chi))), cmap='seismic')

plt.show()
