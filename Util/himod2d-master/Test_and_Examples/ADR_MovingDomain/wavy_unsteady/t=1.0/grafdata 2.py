# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 22:03:54 2012

@author: alonso
@modified: matteo aletti
"""
from numpy import *
from numpy import linalg as LA
import numpy as np
from numpy.linalg.linalg import inv
import matplotlib.pyplot as plt
from pylab import cm, clabel, contour, imshow, colorbar, title, show
from matplotlib.pyplot import pcolor
import os, sys

archi_x  = open('m=4h0.02.out','r')
archi_x1 = open('m=8h0.02.out','r')
archi_x2 = open('m=16h0.02.out','r')
archi_f  = open('ffh002.out','r')

ptosx    = int(archi_x.readline(30))
ptosy    = int(archi_x.readline(30))

ptosx1    = int(archi_x1.readline(30))
ptosy1   = int(archi_x1.readline(30))

ptosx2    = int(archi_x2.readline(30))
ptosy2    = int(archi_x2.readline(30))

ptosxf    = int(archi_f.readline(30))
ptosyf    = int(archi_f.readline(30))

ptos  = ptosx*ptosy
ptos1 = ptosx*ptosy
ptos2 = ptosx*ptosy
ptosf = ptosxf*ptosyf

A = zeros((4,ptos))
A1 = zeros((4,ptos1))
A2 = zeros((4,ptos2))
F = zeros((4,ptosf))

nos = int(sqrt(ptos))

xx = zeros((ptosx,ptosy))
yy = zeros((ptosx,ptosy))
ua = zeros((ptosx,ptosy))

xx1 = zeros((ptosx1,ptosy1))
yy1 = zeros((ptosx1,ptosy1))
ua1 = zeros((ptosx1,ptosy1))

xx2 = zeros((ptosx2,ptosy2))
yy2 = zeros((ptosx2,ptosy2))
ua2 = zeros((ptosx2,ptosy2))

xxf = zeros((ptosxf,ptosyf))
yyf = zeros((ptosxf,ptosyf))
ufem = zeros((ptosxf,ptosyf))

for i in range(ptos):
    archi_x.read(7)
    A[0,i] = archi_x.read(15)
    archi_x.read(5)
    A[1,i] = archi_x.read(15)
    archi_x.read(5)
    A[2,i] = archi_x.readline()

for i in range(ptos1):
    archi_x1.read(7)
    A1[0,i] = archi_x1.read(15)
    archi_x1.read(5)
    A1[1,i] = archi_x1.read(15)
    archi_x1.read(5)
    A1[2,i] = archi_x1.readline()

for i in range(ptos2):
    archi_x2.read(7)
    A2[0,i] = archi_x2.read(15)
    archi_x2.read(5)
    A2[1,i] = archi_x2.read(15)
    archi_x2.read(5)
    A2[2,i] = archi_x2.readline()

for i in range(ptosf):
    archi_f.read(7)
    F[0,i] = archi_f.read(15)
    archi_f.read(5)
    F[1,i] = archi_f.read(15)
    archi_f.read(5)
    F[2,i] = archi_f.readline()


      
for i in range(ptosx):
	xx[i] = A[0,ptosy*i:ptosy+ptosy*i]    
	yy[i] = A[1,ptosy*i:ptosy+ptosy*i]
	ua[i] = A[2,ptosy*i:ptosy+ptosy*i]

for i in range(ptosx1):
	xx1[i] = A1[0,ptosy1*i:ptosy1+ptosy1*i]    
	yy1[i] = A1[1,ptosy1*i:ptosy1+ptosy1*i]
	ua1[i] = A1[2,ptosy1*i:ptosy1+ptosy1*i]
for i in range(ptosx2):
	xx2[i] = A2[0,ptosy2*i:ptosy2+ptosy2*i]    
	yy2[i] = A2[1,ptosy2*i:ptosy2+ptosy2*i]
	ua2[i] = A2[2,ptosy2*i:ptosy2+ptosy2*i]

for i in range(ptosxf):
	xxf[i] = F[0,ptosyf*i:ptosyf+ptosyf*i] 
	yyf[i] = F[1,ptosyf*i:ptosyf+ptosyf*i]
	ufem[i]= F[2,ptosyf*i:ptosyf+ptosyf*i]

archi_x.close()
archi_x1.close()
archi_x2.close()
archi_f.close()

vm=-0.02
vM=0.12
CBt=np.linspace(0.0,vM,16)
CB=np.linspace(vm,vM,16)

plt.subplot(8,1,(1,2))
C=plt.contourf(xx, yy, ua, CB)
#plt.contour(C, levels=C.levels[::1],colors = 'black',hold='on')
plt.ylabel(r'm=4')
plt.yticks(np.arange(-1, 1.4, 0.4))
plt.colorbar(C,use_gridspec=True,ticks=CBt[::3])
plt.grid()

plt.subplot(8,1,(3,4))
C1=plt.contourf(xx1, yy1, ua1, CB)
#plt.contour(C1, levels=C1.levels[::1],colors = 'black',hold='on')
plt.ylabel(r'm=8')
plt.yticks(np.arange(-1, 1.4, 0.4))
plt.colorbar(C,use_gridspec=True,ticks=CBt[::3])
plt.grid()


plt.subplot(8,1,(5,6))
F=plt.contourf(xx2, yy2, ua2, CB)
#plt.contour(F, levels=C.levels[::1],colors = 'black',hold='on')
plt.ylabel(r'm=16')
plt.yticks(np.arange(-1, 1.4, 0.4))
plt.colorbar(C,use_gridspec=True,ticks=CBt[::3])
plt.grid()


plt.subplot(8,1,(7,8))
C2=plt.contourf(xxf, yyf, ufem, CB)
#plt.contour(C2, levels=C.levels[::1],colors = 'black', hold='on')
plt.ylabel(r'FreeFem')

plt.yticks(np.arange(-1, 1.4, 0.4))
plt.colorbar(C,use_gridspec=True,ticks=CBt[::3])
plt.grid()

#plt.show
plt.tight_layout()
plt.savefig('wavy_unsteady.eps',format='eps',transparent='true')
plt.savefig('wavy_unsteadys.png',format='png',transparent='true')
