function [u,p,ta,q,dq,alpha]=PulsatileFlow(r,ru,mu,freq,p0,pn,phi,timestep,grid)
%PulsatileFlow PulsatileFlow function to generate pulsatile flow in tube.
%   u = PulsatileFlow(r,ru,mu,freq,p0,pn,phi,timestep,grid) provides an array u(t,x,y) with t as time step, x and y as the spatial coordinate with size grid x grid.
%   Input parameters for the functions are defined as below:
%
%       r           radius of the tube.
%       ru          density of the fluid in tube.
%       mu          viscosity of the fluid in tube.
%       freq        frequency of oscillation of the fluid.
%       p0          pressure gradient for the steady term of flow. For the pulsatile wave fft this parameter is equal to A0.
%       pn          array of amplitude parameters for the pulsatile wave fft when p(t)=p0+Sum[pn(n)*exp(i*(2*pi*freq*t*n-phi(n))))].
%       phi         array of phase change in the pulsatile wave fft when p(t)=p0+Sum[pn(n)*exp(i*(2*pi*freq*t*n-phi(n))))].
%       timestep    number of time frames to generate the output array. 
%       grid        number of spacial poisitions where data are presented as velocity (size of the array).
%
%
%   phi and pn are arrays of the same size.
%   If either input is a scalar, it is expanded to sinus wave form.
%
%   [u,p,ta,q,dq,alpha] = PulsatileFlow(r,ru,mu,freq,p0,pn,phi,timestep,grid) is the same as u = PulsatileFlow(r,ru,mu,freq,p0,pn,phi,timestep,grid) with more output.
%
%   Output parameters for the functions are defined as below:
%
%       u           u(t,x,y) array with t as time step, x and y as the spatial coordinate including the velocity field.
%       p           p(t) array of pressure data in the domain.
%       ta          ta(t) is an array of the wall shear stress on the wall at different time frames.
%       q           q(t) is an array of the volumetric flow in the tube at different time steps.
%       dq          dq(t) is an array of the volumetric flow rate in the tube at
%       different time steps.
%       alpha       Womersley Number, is a scalar number.
%
%   Examples:
%
%       u = PulsatileFlow(0.001,1060,0.0035,1.4,0,[0.78 1.32 -0.74],[-0.01 -1.45 -0.46],16,32) 
%
%       generates the blood flow in an artery of radius 0.001 m with other
%       parameters set in SI.
%
%       is possible to use SurfMovie(u) function to follow the flow
%       velocity distribution.
%
%   This function uses a Bessel function besselj in Matlab.
%
%   Reference:
%   A. Pashaei and N Fatouraee, "An analytical phantom for the evaluation
%   of medical flow imaging algorithms.", Phys Med Biol. 2009 Mar 21;54(6):1791-821.
%
%   Copyright 2007-2013 Ali Pashaei. 
%   $Revision: 2.0 $  $Date: 2013/11/12 12:31:00 $

ofst=round(grid/2);
rxl=round(3*ofst/4);
h=r/rxl; 
nw=length(pn);
omega=2*pi*freq; u=zeros(timestep,grid,grid); p=zeros(1,timestep); ta=zeros(1,timestep); 
q=zeros(1,timestep); dq=zeros(1,timestep); zt=zeros(1,timestep+1); 
alpha=r*sqrt(omega*ru/mu);
kapa=alpha*i^(1.5)/r;
% to make the waveform
snw=nw*(nw+1)/2; alpha=alpha*sqrt(snw);
for k=1:timestep
    t=k/timestep/freq;
    for l=1:nw
        zt(k)=zt(k)+pn(l)*exp(i*(omega*t*l-phi(l)));
    end
end
% velocity u
CJA = besselj(0,kapa*r);
for m=-rxl:rxl;
    for n=-rxl:rxl;
        for k=1:timestep
            ri=sqrt(m^2+n^2);
            if ri*h<r
                CBJ0=besselj(0,kapa*h*ri);
                u(k,m+ofst,n+ofst)=p0*((ri*h)^2-r^2)/4/mu+real(i/ru/omega/snw*(1-CBJ0/CJA)*zt(k));
            end
        end
    end
end
% pressure, shear rate, volume flow and flow rate
rtst=r*0.93;
CBJ0=besselj(0,kapa*r);
CBJ1=besselj(1,kapa*rtst);
for k=1:timestep
    p(k)=real(zt(k));
    ta(k)=-mu*p0*rtst/2/mu+real(i*kapa/ru/omega/snw*(CBJ1/CBJ0)*zt(k));
    q(k)=p0*pi*(rtst^4/2-r^2*rtst^2)/4/mu+real(i*2*pi/ru/omega/snw*(rtst^2/2-rtst*CBJ1/CBJ0/kapa)*zt(k));
    dq(k)=real(-2*pi/ru*(rtst^2/2-rtst*CBJ1/CBJ0/kapa)*zt(k));
end
end