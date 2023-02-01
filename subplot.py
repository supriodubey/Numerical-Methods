#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:32:27 2022

@author: draco
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

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
           
            i+=1
    f.close()
    return(m1,m2)
##### main ####
mass1, mass2 = readtxt('time_BHillustris1_30.dat')

######## swapping masses #########3
for i in range(len(mass1)):
    if (mass1[i] < mass2[i]):
        mass1[i],mass2[i] =mass2[i],mass1[i]
        
########
N = 100

dm=(max(mass1)-min(mass2))/N

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

##### Chirp Mass #####

mtot = mass1+mass2

for a in range(len(mtot)):
    
    if (mtot[a] <= 0):
        mtot[a]= 0.1
    

mchirp = ((mass2*mass1)**(3/5))/(mtot**(1/5))
print(mchirp)
z2,xedges2,yedges2=np.histogram2d(mtot,mchirp,bins= 50,density=False)
x2=np.zeros(len(xedges2)-1)
y2=np.zeros(len(xedges2)-1)
for i in range(len(xedges2)-1):
    x2[i]=(xedges2[i]+xedges2[i+1])/2
    y2[i]=(yedges2[i]+yedges2[i+1])/2
z2=np.where(z2<=0,z2+1,z2)
z2=np.transpose(z2)

######subplots ########
fig,axs=plt.subplots(1,2)

cs=axs[0].contourf(xedges,yedges,matrix,50,cmap=cm.viridis, norm=LogNorm())
axs[0].set_xlabel('$m_1[M{_{\odot}}$]')
axs[0].set_ylabel('$m_2[M{_{\odot}}$]')


cs2=axs[1].contourf(x2,y2,z2,norm=colors.LogNorm())

plt.tight_layout()
plt.show()
#9223372036314334935 1258065.34091 13.3308464098 0.0 0.0002 742428.0 3.7754 3.134 2.7475 10.5833464098