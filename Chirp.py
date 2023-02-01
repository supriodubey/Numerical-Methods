#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:35:49 2022
excercise - 3
@author: draco
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.colors import LogNorm





chirp_mass =np.genfromtxt("chirpmass_bin.dat",comments="#",usecols=(0), unpack=True)
t_merg = np.genfromtxt("tmerg_bin.dat",comments="#",usecols=(0), unpack=True)
chirp_matrix =np.genfromtxt("chirpmass_tmerg_tot.dat",comments="#")


cs = plt.contourf(t_merg,chirp_mass,chirp_matrix, level = 10., cmap=cm.Reds, norm=LogNorm())

cbar=plt.colorbar(cs,orientation='vertical', label='$N_{merg}$')
cbar.solids.set_edgecolor("face")
cbar.set_label('Number of Merger')

cs=plt.contour(t_merg,chirp_mass,chirp_matrix,[1e4, 5e4, 1e5 ,5e5],linewidths=0.5,  colors='k')
plt.clabel(cs,[1e4, 5e4, 1e5 ,5e5], fontsize=10, inline=True)

plt.xlabel('$t_{merg}$[Gyr]',fontsize=14)
plt.ylabel('$m_{chirp}$[M$_{\odot}$]',fontsize=14)
plt.xlim([0,14])
plt.ylim([0,37])



plt.tight_layout()
plt.show()