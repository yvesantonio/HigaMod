#! /usr/bin/python

from numpy import *
from numpy import linalg as LA
import numpy as np
from numpy.linalg.linalg import inv
import matplotlib.pyplot as plt
from pylab import cm, clabel, contour, imshow, colorbar, title, show
from matplotlib.pyplot import pcolor
import os, sys


def read(filename):
    f=open(filename,'r')
    pointsx = int(f.readline(30))
    pointsy = int(f.readline(30))
    pt = pointsx*pointsy
    A = zeros((4,pt))
    x = zeros((pointsx,pointsy))
    y = zeros((pointsx,pointsy))
    u = zeros((pointsx,pointsy))
    for i in range(pt):
        f.read(7)
        A[0,i] = f.read(15)
        f.read(5)
        A[1,i] = f.read(15)
        f.read(5)
        A[2,i] = f.readline()
    for i in range(pointsx):
	x[i] = A[0,pointsy*i:pointsy+pointsy*i]
        y[i] = A[1,pointsy*i:pointsy+pointsy*i]
        u[i] = A[2,pointsy*i:pointsy+pointsy*i]
    f.close()
    return (x,y,u)
#archi_x = open('himod2.out','r')
#archi_x1 = open('himod4.out','r')
#archi_x2 = open('himod8.out','r')
#archi_f = open('DRFF.out','r')

(xxf,yyf,ufem) = read('DRFF100.out')
(xx,yy,ua)     = read('DRFF200.out')
(xx1,yy1,ua1)  = read('DRFF100.out')
(xx2,yy2,ua2)  = read('DRFF200.out')

vm=0.05
vM=0.13
CB=np.linspace(vm,vM,5)

plt.subplot(8,1,(1,2))
C=plt.contourf(xxf, yyf, ufem, 20,vmin=vm,vmax=vM)
plt.ylabel('h=0.01 (come attualmente nel paper')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()

plt.subplot(8,1,(3,4))
plt.contourf(xx, yy, ua, 20,vmin=vm,vmax=vM)
plt.ylabel('h=0.005')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()


plt.subplot(8,1,(5,6))
plt.contourf(xx1, yy1, ua1, 20,vmin=vm,vmax=vM)
plt.ylabel('m = 4')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()


plt.subplot(8,1,(7,8))
plt.contourf(xx2, yy2, ua2, 20,vmin=vm,vmax=vM)
plt.ylabel('m = 8')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()

plt.show()
#plt.tight_layout()
#plt.savefig('prova.eps',format='eps',transparent='true')
#plt.savefig('prova.png',format='png',transparent='true')
