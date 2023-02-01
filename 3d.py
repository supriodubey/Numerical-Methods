#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:11:21 2022

@author: draco
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors

mettalicity = np.genfromtxt("time_BHillustris1_30.dat",comments="#",usecols=(3), unpack=True,max_rows=10000)
mass1 = np.genfromtxt("time_BHillustris1_30.dat",comments="#",usecols=(6), unpack=True,max_rows=10000)
mass2 = np.genfromtxt("time_BHillustris1_30.dat",comments="#",usecols=(7), unpack=True,max_rows=10000)
tdelay = np.genfromtxt("time_BHillustris1_30.dat",comments="#",usecols=(8), unpack=True,max_rows=10000)

mtot = mass1 + mass2
mchirp = ((mass2*mass1)**(3/5))/(mtot**(1/5))
########## Plot3d #########

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
X = mtot
Y = mettalicity
Z = tdelay
surf=ax.scatter(X,Y,Z,c=mchirp)
ax.set_xlabel('Total Mass (M_sun)')
ax.set_ylabel('Progenitor_s Metallicity')
ax.set_zlabel('Delay time (Gyr)')
fig.colorbar(surf, shrink=0.5, aspect=10)
plt.tight_layout()
plt.show()

