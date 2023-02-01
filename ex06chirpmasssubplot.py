import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

mass1,mass2=np.genfromtxt('time_BHillustris1_30.dat',dtype='float',comments='#',usecols=(6,7),unpack=True,max_rows=10000)
massM=np.where(mass1<mass2,mass2,mass1)
massm=np.where(mass1<mass2,mass1,mass2)

#Preparazione m1-m2
z,xedges,yedges=np.histogram2d(massM,massm,bins=50,density=False)
x=np.zeros(len(xedges)-1)
y=np.zeros(len(xedges)-1)
for i in range(len(xedges)-1):
    x[i]=(xedges[i]+xedges[i+1])/2
    y[i]=(yedges[i]+yedges[i+1])/2
z=np.where(z<=0,z+1,z)
z=np.transpose(z)

#Preparazione mc-mtot
mtot=massM+massm
mc=(massM*massm)**(3/5)/mtot**(1/5)
z2,xedges2,yedges2=np.histogram2d(mtot,mc,bins=50,density=False)
x2=np.zeros(len(xedges2)-1)
y2=np.zeros(len(xedges2)-1)
for i in range(len(xedges2)-1):
    x2[i]=(xedges2[i]+xedges2[i+1])/2
    y2[i]=(yedges2[i]+yedges2[i+1])/2
z2=np.where(z2<=0,z2+1,z2)
z2=np.transpose(z2)


#Subplots
fig,axs=plt.subplots(1,2)

cs=axs[0].contourf(x,y,z,norm=colors.LogNorm())
axs[0].set_xlim(0,45)
axs[0].set_ylim(0.45)

cs2=axs[1].contourf(x2,y2,z2,norm=colors.LogNorm())

plt.tight_layout()
plt.show()
