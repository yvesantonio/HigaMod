#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

dir="./"
def importSol(dir,filename):
    with open(dir+filename,'r') as f:
        L=f.readlines();
        nVisX=int(L[0])
        nVisY=int(L[1])
        x=list()
        y=list()
        s=list()
        for k in range(2,len(L)):
            nbs = L[k].split()
            x.append(float(nbs[0]))
            y.append(float(nbs[1]))
            s.append(float(nbs[2]))
    xx = np.zeros((nVisX,nVisY))
    yy = np.zeros((nVisX,nVisY))
    ss = np.zeros((nVisX,nVisY))
    for i in range(nVisX):
        xx[i]=x[ nVisY*i : nVisY * ( i + 1 ) ];
        yy[i]=y[ nVisY*i : nVisY * ( i + 1 ) ];
        ss[i]=s[ nVisY*i : nVisY * ( i + 1 ) ];
    return (xx,yy,ss)
def plot(name,lev):
    (xr,yr,sr)=importSol(dir,"test."+name)
    plt.contourf(xr,yr,sr,levels=lev)
    
smin = -0.02
smax= 0.1
lev=np.linspace(smin,smax,20)

plt.figure(1)
plt.subplot(121)
plot("ffem",lev)
plt.subplot(122)
plot("himod",lev)
plt.show()
