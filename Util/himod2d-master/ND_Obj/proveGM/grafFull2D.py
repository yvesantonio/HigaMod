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

import os, sys

archi_f = open('./input/Full2D.out','r')

ptosx = int(archi_f.readline(30))
ptosy = int(archi_f.readline(30))

ptos = ptosx*ptosy

F = zeros((4,ptos))
x = zeros((ptosx))
y = zeros((ptosy))
ufem = zeros((ptosx,ptosy))

for i in range(ptos):
   	archi_f.read(7)
   	F[0,i] = archi_f.read(15)
	archi_f.read(5)
	F[1,i] = archi_f.read(15)
	archi_f.read(5)
	F[2,i] = archi_f.readline() 

for i in range(ptosy):
#       print(F[1,i])
         y[i]=F[1,i]
 
for j in range(ptosx):
#       print(F[0,j])
        x[j]=F[0,j*ptosy]

for i in range(ptosx):
#	xx[i]=F[0,ptosy*i:ptosy+ptosy*i]
#	yy[i]=F[1,ptosy*i:ptosy+ptosy*i]
	ufem[i]= F[2,ptosy*i:ptosy+ptosy*i]

archi_f.close()

vm=0.05
vM=0.13
CB=np.linspace(vm,vM,5)
Y,X = np.meshgrid(y,x)
C=plt.contourf(X,Y,ufem, 20,vmin=vm,vmax=vM)
plt.ylabel('FreeFem++')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()
#plt.show
#plt.tight_layout()
plt.savefig('./fotoprova/Full2D.eps',format='eps',transparent='true')
plt.savefig('./fotoprova/Full2D.png',format='png',transparent='true')
