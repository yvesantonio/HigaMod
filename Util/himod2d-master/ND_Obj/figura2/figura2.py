##################################
from numpy import *
from numpy import linalg as LA
import numpy as np
from numpy.linalg.linalg import inv
import matplotlib.pyplot as plt
from pylab import cm, clabel, contour, imshow, colorbar, title, show
import os, sys
#################################
archi_full = open('./input/Full2D.out','r')
ptosx = int(archi_full.readline(30))
ptosy = int(archi_full.readline(30))
ptos = ptosx*ptosy

Full = zeros((4,ptos))
xfull = zeros((ptosx))
yfull = zeros((ptosy))
ufemfull = zeros((ptosx,ptosy))

for i in range(ptos):
   	archi_full.read(7)
   	Full[0,i] = archi_full.read(15)
	archi_full.read(5)
	Full[1,i] = archi_full.read(15)
	archi_full.read(5)
	Full[2,i] = archi_full.readline() 
for i in range(ptosy):
         yfull[i]=Full[1,i]
for j in range(ptosx):
        xfull[j]=Full[0,j*ptosy]
for i in range(ptosx):
	ufemfull[i]= Full[2,ptosy*i:ptosy+ptosy*i]
archi_full.close()
Yfull,Xfull = np.meshgrid(yfull,xfull)
##########################################################################

archi_U2d = open('./input/Krpo2Dgood.out','r')
archi_U1d = open('./input/matricegood.out','r')
archi_x1d = open('./input/xcoordinategood.out','r')
archi_y1d = open('./input/ycoordinategood.out','r')
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
	F[0,i] = archi_U2d.read(13)
	archi_U2d.read(7)
	F[1,i] = archi_U2d.read(13)
	archi_U2d.read(7)
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


# salvo la mesh in direzione y (prima meta')
for i in range(50):
	archi_y1d.read(1)
#	print(archi_y1d.readline())
	y1d[i]=float(archi_y1d.readline())

# salvo la mesh in direzione y (seconda meta')
for i in range(50):
	archi_y1d.read(2)
	y1d[50+i]=float(archi_y1d.readline()) 

# salvo la mesh in y dal 2D
for i in range(ptosy2d):
	y[i]=F[1,i]

# salvo la mesh in x del 2D in entrambi i posti
for j in range(ptosx2d):
	x[j]=F[0,j*ptosy2d]
	x2d[j]=F[0,j*ptosy2d]

# salvo la mesh in x del 1D in entrambi i posti
for i in range(ptosx1d-1):
	archi_x1d.read(1)
	x[ptosx2d+i]=archi_x1d.readline()

archi_U2d.close()
archi_U1d.close()
archi_x1d.close()
archi_y1d.close()

Y,X = np.meshgrid(y,x)


##########################################################################

archi_U2dbad = open('./input/Krpo2Dbad.out','r')
archi_U1dbad = open('./input/matricebad.out','r')
archi_x1dbad = open('./input/xcoordinatebad.out','r')
archi_y1dbad = open('./input/ycoordinatebad.out','r')
# Convenzione: vengono passati il numero di vertici nelle due direzioni
ptosx2dbad = int(archi_U2dbad.readline(30))
ptosy2dbad = int(archi_U2dbad.readline(30))
ptosbad = (ptosx2dbad)*(ptosy2dbad)
Fbad = zeros((4,ptosbad))
# Modifica qui per griglia 1D diversa
ptosx1dbad = 301
ptosy1dbad = 100
F1dbad = zeros((ptosx1dbad,ptosy1dbad))
# Per fare il plot fuso assumiamo di utilizzare lo stesso numero di vertici
# nelle direzione Y. Nella direzione x vengono sovrapposti due vertici.
xbad = zeros(ptosx2dbad+ptosx1dbad-1)
ybad = zeros(ptosy2dbad)

x1dbad = zeros(ptosx1dbad)
y1dbad = zeros(ptosy1dbad)

x2dbad = zeros(ptosx2dbad)
y2dbad = zeros(ptosy2dbad)

ufem1dbad = zeros((ptosx1dbad,ptosy1dbad))
ufem2dbad = zeros((ptosx2dbad,ptosy2dbad))
ufembad   = zeros((ptosx2dbad + ptosx1dbad - 1,ptosy2dbad))

# Salvo temporanemante i valori della soluzione 2d
for i in range(ptosbad):
   	archi_U2dbad.read(7)
	Fbad[0,i] = archi_U2dbad.read(13)
	archi_U2dbad.read(7)
	Fbad[1,i] = archi_U2dbad.read(13)
	archi_U2dbad.read(7)
	Fbad[2,i] = archi_U2dbad.readline() 

# Salvo il pezzo di soluzione 2D nella soluzione totale e in quella solo 2D
for i in range(ptosx2dbad):
	ufembad[i]= Fbad[2,ptosy2dbad*i:ptosy2dbad+ptosy2dbad*i]
	ufem2dbad[i] = Fbad[2,ptosy2dbad*i:ptosy2dbad+ptosy2dbad*i]

# salvo la soluzione 1d sia in ufem che in ufem1d
for j in range(ptosy1dbad):
	archi_U1dbad.read(2)
	for i in range(ptosx1dbad-3):
		ufembad[ptosx2dbad+i,j] = archi_U1dbad.read(6)
		ufem1dbad[i,j] = ufembad[ptosx2dbad+i,j]
		archi_U1dbad.read(3)
	ufembad[ptosx2dbad+ptosx1dbad-2,j]=archi_U1dbad.read(6)
	ufem1dbad[ptosx1dbad-1,j]=ufem[ptosx2dbad+ptosx1dbad-2,j]
	archi_U1dbad.readline()


# salvo la mesh in direzione y (prima meta')
for i in range(50):
	archi_y1dbad.read(1)
#	print(archi_y1dbad.readline())
	y1dbad[i]=float(archi_y1dbad.readline())

# salvo la mesh in direzione y (seconda meta')
for i in range(50):
	archi_y1dbad.read(2)
	y1dbad[50+i]=float(archi_y1dbad.readline()) 

# salvo la mesh in y dal 2D
for i in range(ptosy2dbad):
	ybad[i]=Fbad[1,i]

# salvo la mesh in x del 2D in entrambi i posti
for j in range(ptosx2dbad):
	xbad[j]=Fbad[0,j*ptosy2dbad]
	x2dbad[j]=F[0,j*ptosy2dbad]

# salvo la mesh in x del 1D in entrambi i posti
for i in range(ptosx1dbad-1):
	archi_x1dbad.read(1)
	xbad[ptosx2dbad+i]=archi_x1dbad.readline()

archi_U2dbad.close()
archi_U1dbad.close()
archi_x1dbad.close()
archi_y1dbad.close()

Ybad,Xbad = np.meshgrid(ybad,xbad)
############################################################################
archi_sigma = open('./input/sigma','r')
#archi_conc  = open('./input/conc','r')
xs=zeros(22)
ss=zeros(22)
#cc=zeros(22)
for i in range(22):
	xs[i]=archi_sigma.read(9)
	ss[i]=archi_sigma.readline()
#	archi_conc.read(9)
#	cc[i]=archi_conc.readline()
############################################################################

archi_x  = open('./input/himo52.out','r')
archi_x1 = open('./input/himo52bad.out','r')

ptosx = int(archi_x.readline(30))
ptosy = int(archi_x.readline(30))


ptosxb = int(archi_x1.readline(30))
ptosyb = int(archi_x1.readline(30))

ptos = ptosx*ptosy
ptosb = ptosxb*ptosyb

A = zeros((4,ptos))
A1 = zeros((4,ptosb))

nos = int(sqrt(ptos))

xx = zeros((ptosx,ptosy))
yy = zeros((ptosx,ptosy))
xxb = zeros((ptosxb,ptosyb))
yyb = zeros((ptosxb,ptosyb))
ua = zeros((ptosx,ptosy))
ua1= zeros((ptosx,ptosy))

for i in range(ptos):
    archi_x.read(7)
    A[0,i] = archi_x.read(15)
    archi_x.read(5)
    A[1,i] = archi_x.read(15)
    archi_x.read(5)
    A[2,i] = archi_x.readline()

for i in range(ptosb):
    archi_x1.read(7)
    A1[0,i] = archi_x1.read(15)
    archi_x1.read(5)
    A1[1,i] = archi_x1.read(15)
    archi_x1.read(5)
    A1[2,i] = archi_x1.readline()
      
for i in range(ptosx):
	xx[i] = A[0,ptosy*i:ptosy+ptosy*i]    
	yy[i] = A[1,ptosy*i:ptosy+ptosy*i]
	ua[i] = A[2,ptosy*i:ptosy+ptosy*i]

for i in range(ptosxb):
	xxb[i] = A1[0,ptosyb*i:ptosyb+ptosyb*i]    
	yyb[i] = A1[1,ptosyb*i:ptosyb+ptosyb*i]
	ua1[i] = A1[2,ptosyb*i:ptosyb+ptosyb*i]

archi_x.close()
archi_x1.close()
############################################################################
vm=0.02
vM=0.11


CBt=np.linspace(0.0,vM,16)
CB=np.linspace(vm,vM,16)

plt.subplot(4,4,(1,2))
C2=plt.contourf(Xfull,Yfull,ufemfull, CB)
#plt.ylabel('FreeFem++')
plt.yticks(np.arange(0, 1+0.25, 0.25))
plt.xticks(np.linspace(0,10,6))
#plt.colorbar(C2,use_gridspec=True,ticks=CBt[::3])
plt.grid()

#sigma
plt.subplot(4,4,(3,4))
plt.plot(xs,ss)
plt.plot([4,4],[1.1,2],color='black',lw=1)
plt.plot([2.5,2.5],[1.1,2],color='black',lw=1)
plt.xticks(np.linspace(0,10,6))
plt.yticks(np.linspace(1.1,2,4))
#plt.plot(xs,cc)

#good
plt.subplot(4,4,(5,6))
C1=plt.contourf(X,Y,ufem,CB)
plt.plot([4,4],[0.0,1.0],color='black',lw=1)
#plt.ylabel('1D-reduced model')
plt.yticks(np.arange(0, 1+0.25, 0.25))
plt.xticks(np.linspace(0,10,6))
#plt.colorbar(C1,use_gridspec=True,ticks=CBt[::3])
plt.grid()

#bad
plt.subplot(4,4,(7,8))
C1=plt.contourf(Xbad,Ybad,ufembad,CB)
plt.plot([2.5,2.5],[0.0,1.0],color='black',lw=1)
#plt.ylabel('1D-reduced model')
plt.yticks(np.arange(0, 1+0.25, 0.25))
plt.xticks(np.linspace(0,10,6))
#plt.colorbar(C1,use_gridspec=True,ticks=CBt[::3])
plt.grid()

#HiMod52
#plt.subplot(4,4,(9,10))
#C1=plt.contourf(xx,yy,ua,CB)
#plt.plot([4.,4.],[0.0,1.0],color='black',lw=1)
#plt.ylabel('HiMod, m=5,2')
#plt.colorbar(C1,use_gridspec=True,ticks=CBt[::3])
#plt.grid()

#HiMod52b
#plt.subplot(4,4,(11,12))
#C1=plt.contourf(xxb,yyb,ua1,CB)
#plt.plot([2.5,2.5],[0.0,1.0],color='black',lw=1)
#plt.ylabel('HiMod, m=5,2')
#plt.colorbar(C1,use_gridspec=True,ticks=CBt[::3])
#plt.grid()

plt.tight_layout()

plt.savefig('./fotoprova/Figura2.eps',format='eps',transparent='true')
plt.savefig('./fotoprova/Figura2.png',format='png',transparent='true')

#epstool --copy --bbox Figura2.eps Figura2c.eps