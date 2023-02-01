#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:35:49 2022

@author: draco
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.colors import LogNorm




#m1,m2 =np.genfromtxt("time_BHillustris1_30.dat",comments="#",usecols=(6,7), unpack=True)


def readtxt(fname):
    f=open(fname,"r") #r is read mode
    nline=0
    for i in f: #just count the number of line
        nline+=1
    f=open(fname,"r") #strat again from the fist line
    m1=np.zeros(nline,dtype='float') #inizialize arrays
    m2=np.zeros(nline,dtype='float')
    i=0
    for line in f: #skip comments
        if(line[0]==str("#")):
            continue
        word_list = line.split() #breking the list into line
        if(word_list[0]!=str("#")):
            m1[i]=np.float(word_list[6]) #fill the arrays with masses
            m2[i]=np.float(word_list[7])
            #I tried to order the masses here but the result was strange
            #if(word_list[6]>word_list[7]):
            #    m1[i]=np.float(word_list[6])
            #    m2[i]=np.float(word_list[7])
            #if(word_list[6]<=word_list[7]):
            #    m2[i]=np.float(word_list[6])
            #    m1[i]=np.float(word_list[7])
            i+=1
    f.close()
    return(m1,m2)


#main
mass1, mass2 = readtxt('time_BHillustris1_30.dat')
print(mass1[0],mass2[0])
#print(m1[285],m2[285])
#print(max(m1),max(m1),min(m1),min(m2))

######## swapping masses #########3
for i in range(len(mass1)):
    if (mass1[i] < mass2[i]):
        mass1[i],mass2[i] =mass2[i],mass1[i]
########
#print(len(m1))
N = 50

#xedges = np.linspace(max(m1), min(m1) , N )
#yedges = np.linspace(max(m2), min(m2) , N )
#zedges = np.linspace(1,8,0.875)
#print (zedges)
#print(max(m1),max(m1),min(m1),min(m2))
#print(m1[0],m2[0])
dm=(max(mass1)-min(mass2))/N
#print(dm)
xedges = np.zeros(N)
xedges[0]=0.
dm1 = (np.max(mass1)-np.min(mass1))/N
for j in range (N):
    xedges[j]=xedges[j-1]+dm
   

yedges = np.zeros(N)
yedges[0]=0.
dm2 = (np.max(mass2)-np.min(mass2))/N
for k in range (N):
    yedges[k]=yedges[k-1]+dm

matrix=np.zeros([N,N]) #create

for l in range(len(mass1)): #assigne values
    index1=int(mass1[l]/dm)
    index2=int(mass2[l]/dm)
    matrix[index2-1][index1-1]+=1
    
matrix[np.where(matrix<0.01)]=0.1
#Z = np.transpose(Z)
#print(len(Z))


result=plt.contourf(xedges,yedges,matrix,100,cmap=cm.viridis, norm=LogNorm())
plt.colorbar(result,orientation='vertical',label='$N_{merg}$')
plt.xlabel('$m_1[M{_{\odot}}$]')
plt.ylabel('$m_2[M{_{\odot}}$]')

plt.tight_layout()
plt.show()
