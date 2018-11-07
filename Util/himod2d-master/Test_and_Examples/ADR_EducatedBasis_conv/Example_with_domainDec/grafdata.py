# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 22:03:54 2012

@author: alonso
@modified: matteo aletti
"""
from numpy import *
#from numpy import linalg as LA
#import numpy as np
#from numpy.linalg.linalg import inv
import matplotlib.pyplot as plt
#from pylab import cm, clabel, contour, imshow, colorbar, title, show
#from matplotlib.pyplot import pcolor
#import os, sys

archi_x  = open('ND41.out','r')
archi_x1 = open('RR41.out','r')
archi_x2 = open('RR42.out','r')
archi_x3 = open('ND42.out','r')
archi_f  = open('DDFF.out','r')

ptosx = int(archi_x.readline(30))
ptosy = int(archi_x.readline(30))

archi_x1.readline(30);
archi_x1.readline(30);

archi_x2.readline(30);
archi_x2.readline(30);

archi_x3.readline(30);
archi_x3.readline(30);

archi_f.readline(30)
archi_f.readline(30)

ptos = ptosx*ptosy

A = zeros((4,ptos))
A1 = zeros((4,ptos))
A2 = zeros((4,ptos))
A3 = zeros((4,ptos))
F = zeros((4,ptos))


nos = int(sqrt(ptos))

xx = zeros((ptosx,ptosy))
yy = zeros((ptosx,ptosy))
ua = zeros((ptosx,ptosy))
ua1= zeros((ptosx,ptosy))
ua2= zeros((ptosx,ptosy))
ua3= zeros((ptosx,ptosy))
ufem = zeros((ptosx,ptosy))

for i in range(ptos):
    archi_x.read(7)
    A[0,i] = archi_x.read(15)
    archi_x.read(5)
    A[1,i] = archi_x.read(15)
    archi_x.read(5)
    A[2,i] = archi_x.readline()

for i in range(ptos):
    archi_x1.read(7)
    A1[0,i] = archi_x1.read(15)
    archi_x1.read(5)
    A1[1,i] = archi_x1.read(15)
    archi_x1.read(5)
    A1[2,i] = archi_x1.readline()

for i in range(ptos):
    archi_x2.read(7)
    A2[0,i] = archi_x2.read(15)
    archi_x2.read(5)
    A2[1,i] = archi_x2.read(15)
    archi_x2.read(5)
    A2[2,i] = archi_x2.readline()

for i in range(ptos):
    archi_x3.read(7)
    A3[0,i] = archi_x3.read(15)
    archi_x3.read(5)
    A3[1,i] = archi_x3.read(15)
    archi_x3.read(5)
    A3[2,i] = archi_x3.readline()

for i in range(ptos):
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
	ua1[i]= A1[2,ptosy*i:ptosy+ptosy*i]
	ua2[i]= A2[2,ptosy*i:ptosy+ptosy*i]
	ua3[i]= A3[2,ptosy*i:ptosy+ptosy*i]
	ufem[i]= F[2,ptosy*i:ptosy+ptosy*i]

archi_x.close()
archi_x1.close()
archi_x2.close()
archi_x3.close()
archi_f.close()

vm=0.0
vM=0.055
CB=linspace(vm,vM,20)
CBt=linspace(vm,vM,5)

plt.subplot(5,2,(1,2))
C=plt.contourf(xx, yy, ufem,CB)
plt.title('FreeFem++')
plt.colorbar(C,use_gridspec=True,ticks=CBt)
plt.grid()

plt.subplot(5,2,(3,4))
plt.contourf(xx, yy, ua,CB)
plt.title(r'Neumann-Dirichlet, $m_1=4,m_2=1$')
plt.plot([1.5,1.5],[0,1],linewidth=1,color='black')
plt.colorbar(C,use_gridspec=True,ticks=CBt)
plt.grid()


plt.subplot(5,2,(5,6))
plt.contourf(xx, yy, ua3,CB)
plt.title(r'Neumann-Dirichlet, $m_1=4,m_2=2$')
plt.plot([1.5,1.5],[0,1],linewidth=1,color='black')
plt.colorbar(C,use_gridspec=True,ticks=CBt)
plt.grid()

plt.subplot(5,2,(7,8))
plt.contourf(xx, yy, ua1,CB)
plt.title(r'Robin-Robin, $m_1=4,m_2=1$')
plt.plot([1.5,1.5],[0,1],linewidth=1,color='black')
plt.colorbar(C,use_gridspec=True,ticks=CBt)
plt.grid()


plt.subplot(5,2,(9,10))
plt.contourf(xx, yy, ua2,CB)
plt.title(r'Robin-Robin, $m_1=4,m_2=2$')
plt.plot([1.5,1.5],[0,1],linewidth=1,color='black')
plt.colorbar(C,use_gridspec=True,ticks=CBt)
plt.grid()

#plt.show
plt.tight_layout()
plt.savefig('DD.eps',format='eps',transparent='true')
plt.savefig('DD.png',format='png',transparent='true')
