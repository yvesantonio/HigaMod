# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 22:03:54 2012

@author: davide Baroli
"""

from numpy import *

import numpy as np
import math
import ipdb
import matplotlib.pyplot as plt
from pylab import cm, clabel, contour, imshow, colorbar, title, show
from matplotlib.pyplot import pcolor
from matplotlib import ticker
import os, sys
from matplotlib import cm
import re # now we use a Hammer for string and whitespace!
# Himod
import argparse

parser = argparse.ArgumentParser(description='Comparison results from FF++ and Matlab HiModLab')
parser.add_argument("-inputF","--input-folder", type=str, action='store',  dest='folderI', help="Folder Input") # allows input of the forward read
parser.add_argument("-outputF","--output-folder", type=str, action='store',  dest='folderO', help="Folder Input") # allows input of the forward read
parser.add_argument("-ff","--input-FreeFem++", type=str, action='store',  dest='inputf', help="FreeFem file") # allows input of the forward read
parser.add_argument("-hh","--input-HiModLab", type=str, action='store',  dest='inputh', help="HiModLab file") # allows input of the forward read
parser.add_argument("-hq","--input-QuadratureNode", type=str, action='store',  dest='inputq', help="HiModLab file on Quadrature Node") # allows input of the forward read
parser.add_argument("-N","--number-modes", type=str, action='store',default='7', dest='modes', help="HiModLab file") # allows input of the forward read
parser.add_argument("-o", "--output", type=str, action='store', dest='output', help="Output filename")
args = parser.parse_args()

inputf=args.inputf
inputh=args.inputh
inputq=args.inputq
output=args.output

folderI=args.folderI
folderO=args.folderO

# This part move to parser way

namefolderIn=os.path.join(os.getcwd(),folderI)

namefile=os.path.join(namefolderIn,inputh)

namefile_qn=os.path.join(namefolderIn,inputq)

#namefile=os.path.join(namefolderIn,'himod_m_7.out')
# Freefem output
#namefilef=os.path.join(namefolderIn,'DDFF.out')

namefilef=os.path.join(namefolderIn,inputf)

# Test1 (directory di output)
namefolder=os.path.join(os.getcwd(),folderO)
if not os.path.exists(namefolder):
    try:
        os.mkdir(namefolder)
    except:
        raise Warning("Output Folder exists")

lines = tuple(open(namefile, 'r'))

ptosx=int(lines[0].strip())
ptosy=int(lines[1].strip())

ptos = ptosx*ptosy

A = zeros((3,ptos))


nos = int(sqrt(ptos))

xx = zeros((ptosx,ptosy))
yy = zeros((ptosx,ptosy))
ua = zeros((ptosx,ptosy))

err = zeros((ptosx,ptosy))

for i in range(0,ptos):
    # We jump first 2 lines
    ArrayString=re.split('\s+', lines[i+2].strip())
    A[0,i]=float(ArrayString[0])
    A[1,i]=float(ArrayString[1])
    A[2,i]=float(ArrayString[2])



lines_f = tuple(open(namefilef, 'r'))

#####################
# New for the error #
#####################

ptosx=int(lines[0].strip())
ptosy=int(lines[1].strip())

ptos = ptosx*ptosy

A = np.zeros((5,ptos))


nos = int(math.sqrt(ptos))

xx = np.zeros((ptosx,ptosy))
yy = np.zeros((ptosx,ptosy))
ua = np.zeros((ptosx,ptosy))
uax = np.zeros((ptosx,ptosy))
uay = np.zeros((ptosx,ptosy))

err = np.zeros((ptosx,ptosy))
errx = np.zeros((ptosx,ptosy))
erry = np.zeros((ptosx,ptosy))

lines_qn = tuple(open(namefile_qn, 'r'))
dimQx=int(lines_qn[0].strip())
dimQy=int(lines_qn[1].strip())
print("dimqx="+str(dimQx))

# 720
xq=np.zeros(dimQx)
wxq=np.zeros(dimQx)
mesh_xx=np.zeros(dimQx)
mesh_wx=np.zeros(dimQx)
# 80
mesh_y=np.zeros(dimQy)
mesh_wy=np.zeros(dimQy)

for i in range(0,ptos):
    # We jump first 2 lines
    ArrayString=re.split('\s+', lines[i+2].strip())
    A[0,i]=float(ArrayString[0])
    A[1,i]=float(ArrayString[1])
    A[2,i]=float(ArrayString[2])
    A[3,i]=float(ArrayString[3])
    A[4,i]=float(ArrayString[4])

#Quadrature objects
for i in range(0,dimQx):
    ArrayString=re.split('\s+', lines_qn[i+2].strip())
    xq[i]=float(ArrayString[0])
    wxq[i]=float(ArrayString[1])
    mesh_xx[i]=float(ArrayString[2])
    mesh_wx[i]=float(ArrayString[3])
for i in range(0,dimQy):
    ArrayString=re.split('\s+', lines_qn[i+dimQx+2].strip())
    mesh_y[i]=float(ArrayString[0])
    mesh_wy[i]=float(ArrayString[1])


lines_f = tuple(open(namefilef, 'r'))

#######################################################################

ptosxf=int(lines_f[0].strip())
ptosyf=int(lines_f[1].strip())
ptosf = ptosxf*ptosyf

F = zeros((5,ptosf))
xxf = zeros((ptosxf,ptosyf))
yyf = zeros((ptosxf,ptosyf))


ufem = zeros((ptosxf,ptosyf))
ufemx = np.zeros((ptosxf,ptosyf))
ufemy = np.zeros((ptosxf,ptosyf))

# why ptos << ptosf ?
for i in range(0,ptosf):
    # We jump first 2 lines
    ArrayString=re.split('\s+', lines_f[i+2].strip())
    F[0,i]=float(ArrayString[0])
    F[1,i]=float(ArrayString[1])
    F[2,i]=float(ArrayString[2])
    F[3,i]=float(ArrayString[3])
    F[4,i]=float(ArrayString[4])


for i in range(ptosx):
    xx[i] = A[0,ptosy*i:ptosy+ptosy*i]
    yy[i] = A[1,ptosy*i:ptosy+ptosy*i]
    ua[i] = A[2,ptosy*i:ptosy+ptosy*i]
    uax[i] = A[3,ptosy*i:ptosy+ptosy*i]
    uay[i] = A[4,ptosy*i:ptosy+ptosy*i]

for i in range(ptosxf):
    xxf[i] =F[0,ptosyf*i:ptosyf+ptosyf*i]
    yyf[i] =F[1,ptosyf*i:ptosyf+ptosyf*i]
    ufem[i]=F[2,ptosyf*i:ptosyf+ptosyf*i]
    ufemx[i]=F[3,ptosyf*i:ptosyf+ptosyf*i]
    ufemy[i]=F[4,ptosyf*i:ptosyf+ptosyf*i]

#for idx in np.where(ufem<1e-16)[0]:
#    ufem[idx]=1e-16

#Comparison step grid HiMod and FF
hxf=2*abs(xxf[0][0]-xxf[1][0])
hx=2*abs(xx[0][0]-xx[1][0])

hy= abs(yy[0][2]-yy[0][1])
hyf= abs(yyf[0][2]-yyf[0][1])

print("ratio in mesh_x::%s"%(hx/hxf))
print("ratio in mesh_y :%s"%(hy/hyf))

# mesh in x
xxt=[xx[i][0] for i in range(ptosx)]

xxtf=[xxf[i][0] for i in range(ptosxf)]
yytf=[yyf[0][i] for i in range(ptosyf)]

idxt= (np.abs(xxtf-xx[1][0])).argmin()
idyt= (np.abs(yytf-yy[1][0])).argmin()

# Closest point
assert(hy/hyf ==1, "Same distretization in y-axis")

indexf=[(np.abs(xxtf-xx[i][0])).argmin() for i in range(ptosx)]

# TODO: complete this extension if hy/hyf!=1
#indexff=[((np.abs(xxtf-xx[i][0])).argmin(),np.abs(yytf-yy[0][i]).argmin()) for i in range(ptosx) for j in range(ptosy)]

ufem_c=np.zeros(np.shape(ua))

ufem_cx=np.zeros(np.shape(uax))
ufem_cy=np.zeros(np.shape(uay))

for i in range(ptosx):
    ufemclosest=F[2,ptosyf*indexf[i]:ptosyf+ptosyf*indexf[i]]
    err[i]=abs(ua[i]-ufemclosest)
    ufem_c[i,:]=ufemclosest
    
    ufemxclosest=F[3,ptosyf*indexf[i]:ptosyf+ptosyf*indexf[i]]
    errx[i]=abs(uax[i]-ufemxclosest)
    ufem_cx[i,:]=ufemxclosest
    
    ufemyclosest=F[4,ptosyf*indexf[i]:ptosyf+ptosyf*indexf[i]]
    erry[i]=abs(uay[i]-ufemyclosest)
    ufem_cy[i,:]=ufemyclosest


############################
# Computation of the error #
############################

ERR    = err*err
ERRx   = errx*errx
ERRy   = erry*erry
ERRtot = ERR+ERRx+ERRy

UFEM  = ufem_c * ufem_c

UFEMx = ufem_cx * ufem_cx

UFEMy = ufem_cy * ufem_cy

UFEMtot = UFEM + UFEMx + UFEMy

#print("mesh_wx :%s"%mesh_wx)
print("ptosx :%s"%ptosx)
print("ptosy :%s"%ptosy)


L2err = math.sqrt(np.inner(np.matmul(ERR.T,mesh_wx),mesh_wy))
H1err =  math.sqrt( np.inner(np.matmul(ERRtot.T,mesh_wx),mesh_wy))
#L2errnorm=L2err/math.sqrt(np.inner(np.matmul(UFEM.T,mesh_wx),mesh_wy))
#H1errnorm=H1err/math.sqrt( np.inner(np.matmul(UFEMtot.T,mesh_wx),mesh_wy))

print("L2 :%e"%L2err)
print("H1 :%e"%H1err)
#print("L2norm :%e"%L2errnorm)
#print("H1norm :%e"%H1errnorm)

#########################
# CREATION OF THE PLOTS #
#########################

# LIMITS OF THE COLOR BAR FOR THE PLOTS

v_min = np.min(ua) * 1.05;
v_max = np.max(ufem) * 1.05;

# CREATE THE PLOTS

CS_f = plt.contourf(xxf, yyf, ufem, 100, cmap=cm.jet,levels=np.linspace(v_min,v_max,50))

CS_h = plt.contourf(xx, yy, ua, 100, cmap=cm.jet,levels=np.linspace(v_min,v_max,50))

CS_e = plt.contourf(xx, yy, err, 100, cmap=cm.jet,levels=np.linspace(v_min,v_max,50))

CS_hx = plt.contourf(xx, yy, uax, 100, cmap=cm.jet,levels=np.linspace(v_min,v_max,50))

# PRINT THE IMAGES OF THE PLOTS

plt.subplot(4,1,(1,2))
plt.contourf(xxf, yyf, ufem, 50,cmap=cm.jet,levels=np.linspace(v_min,v_max,50))
plt.ylabel('FreeFem++')
plt.colorbar(CS_f)
plt.grid()


plt.subplot(4,1,(3,4))
plt.contourf(xx, yy, ua, 50,cmap=cm.jet,levels=np.linspace(v_min,v_max,50))
plt.ylabel('HigaMod')
plt.colorbar(CS_h)
plt.grid()

"""
plt.subplot(8,1,(5,6))
plt.contourf(xx, yy, err )
plt.ylabel('error')
cbar = plt.colorbar(CS_e)
plt.grid()
"""


"""
plt.subplot(8,1,(5,6))
plt.contourf(xx, yy, ua1, 20,vmin=vm,vmax=vM)
plt.ylabel('m = 4')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()


plt.subplot(8,1,(7,8))
plt.contourf(xx, yy, ua2, 20,vmin=vm,vmax=vM)
plt.ylabel('m = 8')
plt.colorbar(C,use_gridspec=True,ticks=CB)
plt.grid()
"""
plt.tight_layout()
plt.savefig(namefolder+'/'+output+'.eps',format='eps',transparent='true')
plt.savefig(namefolder+'/'+output+'.pdf',format='pdf',transparent='false') #'true'
#plt.show()
#plt.savefig(namefolder+'/DD.eps',format='eps',transparent='true')
#plt.savefig(namefolder+'/DD.png',format='png',transparent='true')

############################
# Write file with the error
############################

L2err = round(L2err,6);
H1err = round(H1err,6);

outputR=np.asarray(['L2='+str(L2err),'H1='+str(H1err),'hx='+str(hx), 'm='+args.modes])
print(outputR)

if os.path.exists(namefolder+'/'+output+'.txt'):
    with open(namefolder+'/'+output+'.txt', 'a') as fileout:
        for item in outputR:
            fileout.write("%s\n" % item)
else:
    
    with open(namefolder+'/'+output+'.txt', 'wb') as fileout:
        for item in outputR:
            fileout.write("%s\n" % item)
