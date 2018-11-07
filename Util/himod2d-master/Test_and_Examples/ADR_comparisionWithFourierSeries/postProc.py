#! /usr/bin/python

import os.path as pt
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

parser = argparse.ArgumentParser(description='This is a post processing tool')
parser.add_argument('files', type=str, nargs='+', help='file to plot')
parser.add_argument('-t','--titles', type=str, nargs='*', help='titles')
parser.add_argument('-s','--subplot',type=int, nargs=2,default=[1, 1], help='subplot code', required=False)
parser.add_argument('-r','--crange',type=float, nargs=2,default=[0, 1], help='max color range for contourf', required=False)
parser.add_argument('-n','--nclines',type=int,default=20, help='number of contourlines', required=False)
parser.add_argument('-v','--verbosity',type=int,default=0, help='the verbosity',required=False)
args = parser.parse_args()

dir="./"
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
    (x,y,s)=importSol(dir,name)
    if args.verbosity > 0:
        M=s.max()
        m=s.min()
        print name + " min = " + str(m) + " max = " + str(M)
    plt.contourf(x,y,s,levels=lev)

smin = args.crange[0]
smax= args.crange[1]
lev=np.linspace(smin,smax,args.nclines)

plt.figure(1)
c=1;
if args.subplot[0]*args.subplot[1] < len(args.files):
    args.subplot[0]=len(args.files);
    args.subplot[1]=1;
    
cLocToRow=1;
for file in args.files:
	subplot=args.subplot[0]*100+args.subplot[1]*10+c;
	plt.subplot(subplot, aspect='equal')
	plot(file,lev)
        if args.titles is not None and len(args.titles) >= c:
            plt.title(args.titles[c-1])
        else:
            plt.title(file)
        if cLocToRow == args.subplot[1]:
            divider = make_axes_locatable(plt.gca())
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(cax=cax)
            cLocToRow=0
	c=c+1
        cLocToRow=cLocToRow+1
plt.show()
