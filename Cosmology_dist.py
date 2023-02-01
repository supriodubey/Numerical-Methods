#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:18:51 2022

@author: draco
excercise 1 and 2.
"""
import numpy as np
import scipy.integrate as integrate
from astropy.cosmology import WMAP9 as cosmo
from astropy import constants as const
import matplotlib.pyplot as plt

class constant:
    omegaM = 0.2726
    omegaL = 0.7274

def integrand(x,disttype):
    f = 0.
    if (disttype =="comoving"):
        f = 1./(constant.omegaM*(1.+x)**3.+ constant.omegaL)**0.5
    elif (disttype == "lookback" ):
       f=1./((1.+x)*(constant.omegaM*(1.+x)**3.+ constant.omegaL)**0.5)
    else:
        print ("recheck integrand")
    
    return f

def comovingdist(z):
    D_c = integrate.quad(integrand, 0.0, z,"comoving")
    D_comoving = D_c[0]
    D_comoving *= const.c.to('km/s')/cosmo.H(0)         
    
    return D_comoving

def luminositydist(z):
    
    D_lumin = (1.+z)*comovingdist(z)
    return D_lumin

def lookbacktime(z):
    L_t = integrate.quad(integrand, 0.0, z,"lookback")
    T_lookback = L_t[0]
    T_lookback /= cosmo.H(0)
    
    return T_lookback

def properdist(z):
    
    D_proper = const.c.to('km/s')*lookbacktime(z)
    return D_proper

########## MAIN ##########

z = 9.
redshift = []
cd = []
lt = []
ld = []
pd = []

while (z > 0.0):
    redshift.append(z)
    cd.append(comovingdist(z))
    lt.append(lookbacktime(z))
    ld.append(luminositydist(z))
    pd.append(properdist(z))
    
    z -= 0.1
    
#### writing the files #######

fname = "distance.txt"
f =open(fname,"w") 
f.write ("# redshift,Comoving,Lookback,Luminosity,Proper \n") 

for i in range (len(redshift)):
    if (redshift[i]!= 0.0):
         f.write(str(redshift[i])+" "+str(cd[i])+" "+str(lt[i])+" "+str(ld[i])+" "+str(pd[i])+"\n")
###### read ####
fname = "distance.txt"
f = open(fname, "r")
redshift,cd,lt,ld,pd= np.genfromtxt(fname,dtype="float", \
comments="#", usecols=(0,1,3,8,10), unpack=True)

######## plot #####
plt.xlabel("Redshift $z$")
plt.ylabel("Distance{Mpc}")
plt.plot(redshift,cd,color = 'red', linestyle = '--', label = "Comvoving Distance")
plt.plot(redshift,ld, label="luminosity distance")
plt.plot(redshift,pd, color = 'green', linestyle = '-.', label = "Proper Distance")
plt.plot(redshift,lt, color = 'black', linestyle = ':',label = "Lookback time")
plt.legend()
plt.show()
  