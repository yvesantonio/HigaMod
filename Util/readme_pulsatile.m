
% Toolkit for Pulsatile flow in arteries
%
% IMPORTANT > This toolkit is using some functions of Image processing 
% toolkit in matlab. Image processing toolkit is used to export dicom files 
% generated from the toolkit. If you need to export dicom files 
% this toolkit should also be installed.
%
% This toolkit includes a list of functions to provide pulsatile flow images  
% acquired from Womersley solution in the rigid tube. The toolkit includes  
% functions for exact solution of the domain for flow characteristics such
% as velocity, pressure, wall shear rate, volumetric flow and flow rate in
% the computational domain. Also there is a function to generate PC-MR
% images as numerical phantom for the domain. Fourier series for generation
% of temporal wave form are supported when pressure is given as: 
% p(t)=p0+Sum[pn(n)*exp(i*(2*pi*freq*t*n-phi(n))))].
%
%
%******************** List of functions 
%
% [u,p,dudr,q,dq,alpha]=PulsatileFlow(r,ru,mu,freq,p0,pn,phi,timestep,grid);
%
% [prefix,num,venc,pixel_size] = PulsatileFlowDICOM(r,ru,mu,freq,p0,pn,phi,timestep,grid,folder);
%
%
%   Copyright 2007-2013 Ali Pashaei. 
%   $Revision: 3.0 $  $Date: 2013/11/12 11:22:10 $

clf; clear all; clc;

r=0.001; grid=64; timestep=28; 
ru=1060; freq=1.0; mu=0.0035; 

%p0=0.0; pn=[7.99 4.43 0.99 0.74]; phi=[3.1219 -1.6675 0.4821 0.3938]; 
p=[-7.7183,-8.2383,-8.6444,-8.8797,-9.6337,-10.5957,-11.8705,-10.0942,-6.2839,-1.1857,2.6043,4.4323,6.1785,7.8211,9.1311,9.9138,10.3447,10.4011,10.2807,9.8951,8.0597,5.6717,2.5232,1.3301,1.4405,1.9094,1.8145,0.8738,0.7055,0.7343,0.7788,0.7495,0.6711,-0.4796,-1.6541,-2.8643,-3.4902,-4.1714,-5.6581,-6.8024];
[p0,pn,phi] = FourierSeries(p,0.025,4,1/freq);

[u,p,ta,q,dq,alpha]=PulsatileFlow(r,ru,mu,freq,p0,pn,phi,timestep,grid);
plot(q)
SurfMovie(u);
plot (p)
folder='test1';
[prefix,num,venc,pixel_size] = FlowDICOM(u,r,folder)
