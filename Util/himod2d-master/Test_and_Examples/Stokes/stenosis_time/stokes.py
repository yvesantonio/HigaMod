from numpy import *
from numpy import linalg as LA
import numpy as np
from numpy.linalg.linalg import inv
import matplotlib.pyplot as plt
from pylab import cm, clabel, contour, imshow, colorbar, title, show
from matplotlib.pyplot import pcolor
import os, sys

name='m16N20t2'
archi_x = open(name+'.out','r')

ptosx = int(archi_x.readline(30))
ptosy = int(archi_x.readline(30))

ptos = ptosx*ptosy

UXm = zeros((4,ptos))
UYm = zeros((4,ptos))
Pm  = zeros((4,ptos))

xx = zeros((ptosx,ptosy))
yy = zeros((ptosx,ptosy))
ux = zeros((ptosx,ptosy))
uy = zeros((ptosx,ptosy))
p  = zeros((ptosx,ptosy))

for i in range(ptos):
    archi_x.read(7)
    UXm[0,i] = archi_x.read(15)
    archi_x.read(5)
    UXm[1,i] = archi_x.read(15)
    archi_x.read(5)
    UXm[2,i] = archi_x.readline()

for i in range(ptos):
    archi_x.read(7)
    UYm[0,i] = archi_x.read(15)
    archi_x.read(5)
    UYm[1,i] = archi_x.read(15)
    archi_x.read(5)
    UYm[2,i] = archi_x.readline()

for i in range(ptos):
    archi_x.read(7)
    Pm[0,i] = archi_x.read(15)
    archi_x.read(5)
    Pm[1,i] = archi_x.read(15)
    archi_x.read(5)
    Pm[2,i] = archi_x.readline()
      
for i in range(ptosx):
	xx[i] = UXm[0,ptosy*i:ptosy+ptosy*i]    
	yy[i] = UXm[1,ptosy*i:ptosy+ptosy*i]
	ux[i] = UXm[2,ptosy*i:ptosy+ptosy*i]
	uy[i] = UYm[2,ptosy*i:ptosy+ptosy*i]
	p[i]  = Pm[2,ptosy*i:ptosy+ptosy*i]
archi_x.close()

Yt=np.linspace(-1,1,5);

CUX=np.linspace(-0.0015,0.0002,15);
CUXt=np.linspace(-0.0015,0.0002,11);
plt.subplot(3,2,(1,2))
CX=plt.contourf(xx, yy, ux,CUX)
plt.ylabel(r'$U_x$')
plt.yticks(Yt)
plt.colorbar(CX,use_gridspec=True,ticks=CUXt)
plt.grid()

CUY=np.linspace(-0.0007,0.0007,10);
CUYt=np.linspace(-0.0007,0.0007,11);
plt.subplot(3,2,(3,4))
CY=plt.contourf(xx, yy, uy,CUY)
plt.ylabel(r'$U_y$')
plt.yticks(Yt)
plt.colorbar(CY,use_gridspec=True,ticks=CUYt)
plt.grid()

CP=np.linspace(-0.0002,0.0002,15)
CPt=np.linspace(-0.0002,0.0002,11)
plt.subplot(3,2,(5,6))
P=plt.contourf(xx, yy, p,CP)
plt.ylabel('P')
plt.yticks(Yt)
plt.colorbar(P,use_gridspec=True,ticks=CPt)
plt.grid()

#plt.show
plt.tight_layout()
plt.savefig(name+'.eps',format='eps',transparent='true')
plt.savefig('stokes.png',format='png',transparent='true')
