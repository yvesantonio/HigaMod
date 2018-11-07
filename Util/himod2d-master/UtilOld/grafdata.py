# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 22:03:54 2012

@author: alonso
"""
from numpy import *
from numpy import linalg as LA
import numpy as np
from numpy.linalg.linalg import inv
import matplotlib.pyplot as plt
from pylab import cm, clabel, contour, imshow, colorbar, title, show
from matplotlib.pyplot import pcolor
import os, sys

archi_x = open('/home/matteo/Documenti/HIMOD/CODICI/HIMOD/Hi-Mod2D/Test_and_Examples/Testcases/cfr_polynomial_adv/1/cfr_quadrature_rule/m4b4h001.out','r')
archi_x1 = open('/home/matteo/Documenti/HIMOD/CODICI/HIMOD/Hi-Mod2D/Test_and_Examples/Testcases/cfr_polynomial_adv/1/cfr_quadrature_rule/m8b4h001.out','r')
archi_x2 = open('/home/matteo/Documenti/HIMOD/CODICI/HIMOD/Hi-Mod2D/Test_and_Examples/Testcases/cfr_polynomial_adv/1/cfr_quadrature_rule/m16b4h001.out','r')
archi_x3 = open('/home/matteo/Documenti/HIMOD/CODICI/HIMOD/Hi-Mod2D/Test_and_Examples/Testcases/cfr_polynomial_adv/1/cfr_quadrature_rule/m32b4h001.out','r')
#archi_f = open('/home/matteo/Documenti/HIMOD/CODICI/HIMOD/Hi-Mod2D/Test_and_Examples/Testcases/cfr_polynomial/u_fem.out','r')

ptosx = int(archi_x.readline(30))
ptosy = int(archi_x.readline(30))

#archi_f.readline(30)
#archi_f.readline(30)
archi_x1.readline(30);
archi_x1.readline(30);

archi_x2.readline(30);
archi_x2.readline(30);

archi_x3.readline(30);
archi_x3.readline(30);

ptos = ptosx*ptosy

A = zeros((4,ptos))
A1 = zeros((4,ptos))
A2 = zeros((4,ptos))
A3 = zeros((4,ptos))


nos = int(sqrt(ptos))

xx = zeros((ptosx,ptosy))
yy = zeros((ptosx,ptosy))
ua = zeros((ptosx,ptosy))
ua1= zeros((ptosx,ptosy))
ua2= zeros((ptosx,ptosy))
ua3= zeros((ptosx,ptosy))

#ufem = zeros((ptosx,ptosy))

for i in range(ptos):
    archi_x.read(7)
    A[0,i] = archi_x.read(15)
#    print A[0,i]
    archi_x.read(5)
    A[1,i] = archi_x.read(15)
#    print A[1,i]
    archi_x.read(5)
    A[2,i] = archi_x.readline()
#    print A[2,i]
#    archi_x.readline()
    
#    archi_f.readline()
#    archi_f.readline()
#    A[3,i] = archi_f.readline()

for i in range(ptos):
    archi_x1.read(7)
    A1[0,i] = archi_x1.read(15)
#    print A[0,i]
    archi_x1.read(5)
    A1[1,i] = archi_x1.read(15)
#    print A[1,i]
    archi_x1.read(5)
    A1[2,i] = archi_x1.readline()

for i in range(ptos):
    archi_x2.read(7)
    A2[0,i] = archi_x2.read(15)
#    print A[0,i]
    archi_x2.read(5)
    A2[1,i] = archi_x2.read(15)
#    print A[1,i]
    archi_x2.read(5)
    A2[2,i] = archi_x2.readline()
    

for i in range(ptos):
    archi_x3.read(7)
    A3[0,i] = archi_x3.read(15)
#    print A[0,i]
    archi_x3.read(5)
    A3[1,i] = archi_x3.read(15)
#    print A[1,i]
    archi_x3.read(5)
    A3[2,i] = archi_x3.readline()
      
for i in range(ptosx):
	xx[i] = A[0,ptosy*i:ptosy+ptosy*i]    
	yy[i] = A[1,ptosy*i:ptosy+ptosy*i]
	ua[i] = A[2,ptosy*i:ptosy+ptosy*i]
	ua1[i]= A1[2,ptosy*i:ptosy+ptosy*i]
	ua2[i]= A2[2,ptosy*i:ptosy+ptosy*i]
	ua3[i]= A3[2,ptosy*i:ptosy+ptosy*i]

#	ufem[i] = A[3,ptosy*i:ptosy+ptosy*i]
	
archi_x.close()
archi_x1.close()
archi_x2.close()
archi_x3.close()

#archi_f.close()

plt.subplot(8,1,(1,2))
plt.contour(xx, yy, ua, 20)
p = plt.pcolor(xx,yy,ua,vmin=-0.000004,vmax=0.000025)
plt.ylabel('m = 4')
plt.colorbar()
plt.grid()


plt.subplot(8,1,(3,4))
plt.contour(xx, yy, ua1, 100)
p = plt.pcolor(xx,yy,ua1,vmin=-0.000004,vmax=0.000025)
plt.ylabel('m = 8')
plt.colorbar()
plt.grid()


plt.subplot(8,1,(5,6))
plt.contour(xx, yy, ua2, 20)
p = plt.pcolor(xx,yy,ua2,vmin=-0.000004,vmax=0.000025)
plt.ylabel('m = 16')
plt.colorbar()
plt.grid()

 
plt.subplot(8,1,(7,8))
plt.contour(xx, yy, ua3, 20)
p = plt.pcolor(xx,yy,ua3,vmin=-0.000004,vmax=0.000025)
plt.ylabel('m = 32')
plt.colorbar()
plt.grid()

#plt.subplot(6,1,(3,4))
#plt.contour(xx, yy, ua, 100)
#plt.pcolor(xx,yy,ufem)
#plt.ylabel('Solucao aprox com FreeFem')
#plt.colorbar()
#plt.grid()

#plt.subplot(6,1,(5,6))
##plt.contour(xx, yy, ua, 100)
#plt.pcolor(xx,yy,abs(ua - ufem))
#plt.ylabel('Erro')
#plt.colorbar()
#plt.grid()

plt.show()
