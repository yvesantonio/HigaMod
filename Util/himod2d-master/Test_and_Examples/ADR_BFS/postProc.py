#! /usr/bin/python

import os.path as pt
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
parser = argparse.ArgumentParser(description='This is a post processing tool')
parser.add_argument('-d','--dir', help='directory where the files are located',required=True)
parser.add_argument('-r','--crange',type=float, nargs=2,default=[0, 1], help='max color range for contourf', required=False)
parser.add_argument('-n','--nclines',type=int,default=20, help='number of contourlines', required=False)
parser.add_argument('-v','--verbosity',type=int,default=0, help='the verbosity',required=False)
args = parser.parse_args()

dir="./"+args.dir+"/"
def importSol(dir,filename):
    if pt.isfile(dir+filename):
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
    else:
        xx=np.zeros((2,2));
        yy=xx;
        ss=xx;
    return (xx,yy,ss)
def plot(name,lev):
    (xr,yr,sr)=importSol(dir,"right."+name)
    (xl,yl,sl)=importSol(dir,"left."+name)
    if args.verbosity > 0:
        M=max(sr.max(),sl.max())
        m=min(sr.min(),sl.min())
        print name + " min = " + str(m) + " max = " + str(M)
    plt.contourf(xl,yl,sl,levels=lev)
    plt.contourf(xr,yr,sr,levels=lev)

smin = args.crange[0]
smax= args.crange[1]
lev=np.linspace(smin,smax,args.nclines)

plt.figure(1)
plt.subplot(121)
plot("ffem",lev)
plt.subplot(122)
plot("himod",lev)
plt.colorbar()
plt.show()
