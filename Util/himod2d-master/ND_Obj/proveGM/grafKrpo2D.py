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

archi_U2d = open('./input/Krpo2D.out','r')
archi_U1d = open('./input/matrice.out','r')
archi_x1d = open('./input/xcoordinate.out','r')
archi_y1d = open('./input/ycoordinate.out','r')

# Convenzione: vengono passati il numero di vertici nelle due direzioni
ptosx2d = int(archi_U2d.readline(30))
ptosy2d = int(archi_U2d.readline(30))
ptos = (ptosx2d)*(ptosy2d)
F = zeros((4,ptos))

# Modifica qui per griglia 1D diversa
ptosx1d = 301
ptosy1d = 100
F1d = zeros((ptosx1d,ptosy1d))

# Per fare il plot fuso assumiamo di utilizzare lo stesso numero di vertici
# nelle direzione Y. Nella direzione x vengono sovrapposti due vertici.
x = zeros(ptosx2d+ptosx1d-1)
y = zeros(ptosy2d)

x1d = zeros(ptosx1d)
y1d = zeros(ptosy1d)

x2d = zeros(ptosx2d)
y2d = zeros(ptosy2d)

ufem1d = zeros((ptosx1d,ptosy1d))
ufem2d = zeros((ptosx2d,ptosy2d))
ufem   = zeros((ptosx2d + ptosx1d - 1,ptosy2d))

# Salvo temporanemante i valori della soluzione 2d
for i in range(ptos):
   	archi_U2d.read(7)
#	print(archi_f.read(13))
	F[0,i] = archi_U2d.read(13)
	archi_U2d.read(7)
#	print(archi_f.read(13))
	F[1,i] = archi_U2d.read(13)
	archi_U2d.read(7)
#	print(archi_f.readline())
	F[2,i] = archi_U2d.readline() 

# Salvo il pezzo di soluzione 2D nella soluzione totale e in quella solo 2D
for i in range(ptosx2d):
	ufem[i]= F[2,ptosy2d*i:ptosy2d+ptosy2d*i]
	ufem2d[i] = F[2,ptosy2d*i:ptosy2d+ptosy2d*i]

# salvo la soluzione 1d sia in ufem che in ufem1d
for j in range(ptosy1d):
	archi_U1d.read(2)
	for i in range(ptosx1d-3):
		ufem[ptosx2d+i,j] = archi_U1d.read(6)
		ufem1d[i,j] = ufem[ptosx2d+i,j]
		archi_U1d.read(3)
	ufem[ptosx2d+ptosx1d-2,j]=archi_U1d.read(6)
	ufem1d[ptosx1d-1,j]=ufem[ptosx2d+ptosx1d-2,j]
	archi_U1d.readline()


# salvo la mesh in direzione y (prima metà)
for i in range(50):
	archi_y1d.read(1)
#	print(archi_y1d.readline())
	y1d[i]=float(archi_y1d.readline())

# salvo la mesh in direzione y (seconda metà)
for i in range(50):
	archi_y1d.read(2)
#	print(archi_y1d.readline())
	y1d[50+i]=float(archi_y1d.readline()) 

# salvo la mesh in y dal 2D
for i in range(ptosy2d):
#	print(F[1,i])
	y[i]=F[1,i]

# salvo la mesh in x del 2D in entrambi i posti
for j in range(ptosx2d):
#	print(F[0,j])
	x[j]=F[0,j*ptosy2d]
	x2d[j]=F[0,j*ptosy2d]

# salvo la mesh in x del 1D in entrambi i posti
for i in range(ptosx1d-1):
	archi_x1d.read(2)
#	print(archi_x1d.readline())
	x[ptosx2d+i]=archi_x1d.readline()


archi_U2d.close()
archi_U1d.close()

vm=0.05
vM=0.13	
CB=np.linspace(vm,vM,5)


#Figura completa
Y,X = np.meshgrid(y,x)
C=plt.contourf(X,Y,ufem, 20,vmin=vm,vmax=vM)
plt.ylabel('Krpo')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()
#plt.show
#plt.tight_layout()
plt.savefig('./fotoprova/Krpo.eps',format='eps',transparent='true')
plt.savefig('./fotoprova/Krpo.png',format='png',transparent='true')

#Figura solo 2D
Y2d,X2d = np.meshgrid(y2d,x2d)
C2d=plt.contourf(X2d,Y2d,ufem2d,20,vmin=vm,vmax=vM)
plt.ylabel('Krpo only 2D side')
plt.colorbar(C2d,use_gridspec=True,ticks=CB)
plt.grid()
#plt.show
#plt.tight_layout()
plt.savefig('./fotoprova/Krpo2Donly.eps',format='eps',transparent='true')
plt.savefig('./fotoprova/Krpo2Donly.png',format='png',transparent='true')
