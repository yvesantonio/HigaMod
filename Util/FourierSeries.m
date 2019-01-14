function [p0,pn,phi] = FourierSeries(p,dt,N,T)
%FourierSeries FourierSeries function to generate coeficients of Fourier series for a given discrete data. 
%   [p0,pn,phi] = FourierSeries(p,dt,N,T) provides a constant p0 and two
%   arrays pn and phi with length N within period of T. N is the number of
%   coefficients for the Fourier serie.
%   Input parameters for the functions are defined as below:
%
%       p           array of discrete amplitudes for the pulsatile wave.
%       dt          time step size in which p data are presented.
%       N           Number of coefficients required from the function.
%       T           period of oscillation of the fluid.
%
%
%   p is an array, dt is a double number, N is integer and N is double.
%
%   given array p, [p0,pn,phi] = FourierSeries(p,dt,N,T) provides the
%   coefficients for function p(t)=p0+Sum[pn(n)*exp(i*(2*pi*t*n/T-phi(n))))]
%
%   Output parameters for the functions are defined as below:
%
%       p0          a value showing the steady-state component of peiodic wave.
%       pn          an array of size N presenting the absolute values of Fourier serie.
%       phi         an array of size N presenting the phase angle in Fourier serie.
%
%   Examples:
%
%       p=[-7.7183,-8.2383,-8.6444,-8.8797,-9.6337,-10.5957,-11.8705,-10.0942,-6.2839,-1.1857,2.6043,4.4323,6.1785,7.8211,9.1311,9.9138,10.3447,10.4011,10.2807,9.8951,8.0597,5.6717,2.5232,1.3301,1.4405,1.9094,1.8145,0.8738,0.7055,0.7343,0.7788,0.7495,0.6711,-0.4796,-1.6541,-2.8643,-3.4902,-4.1714,-5.6581,-6.8024];
%       [p0,pn,phi] = FourierSeries(p,0.025,4,1)
%
%       generates the Fourier coefficients for wave p, presenting data in intervals of 0.025 within a period of 1 sec.
%       parameters set in SI.
%
%       is possible to visualize the results given function below:
%
%           tt=0:dt:T;
%           pt=p0;
%           for n=1:N
%               pt=pt+pn(n)*cos(2*pi*tt*n/T-phi(n));
%           end
%           plot (tt,pt,'-k',0:dt:T-dt,p,'x')
%
%
%   Reference:
%   A. Pashaei and N Fatouraee, "An analytical phantom for the evaluation
%   of medical flow imaging algorithms.", Phys Med Biol. 2009 Mar 21;54(6):1791-821.
%
%   Copyright 2007-2013 Ali Pashaei. 
%   $Revision: 3.0 $  $Date: 2013/11/12 12:31:00 $

if nargin == 0
    p=[-7.7183,-8.2383,-8.6444,-8.8797,-9.6337,-10.5957,-11.8705,-10.0942,-6.2839,-1.1857,2.6043,4.4323,6.1785,7.8211,9.1311,9.9138,10.3447,10.4011,10.2807,9.8951,8.0597,5.6717,2.5232,1.3301,1.4405,1.9094,1.8145,0.8738,0.7055,0.7343,0.7788,0.7495,0.6711,-0.4796,-1.6541,-2.8643,-3.4902,-4.1714,-5.6581,-6.8024];
    dt = 0.025;
    N = 10;
    T = 1;
    [p0,pn,phi] = FourierCoefficients2(p,dt,N,T);

    tt=0:0.001:1;
    pt=p0;
    for n=1:N
        pt=pt+pn(n)*cos(2*pi*tt*n/T-phi(n));
    end
    plot (tt,pt,'-k',0:dt:T-dt,p,'x')
end
[p0,pn,phi] = FourierCoefficients2(p,dt,N,T);

end

function [A0,AN,BN] = FourierCoefficients(p,dt,N,T)
t=0:dt:T-dt;
A0=0;
for n=1:length(p)
    A0=A0+p(n)*dt/T;
end
AN=zeros(N,1);
BN=zeros(N,1);
for n=1:N;
    for m=1:length(p)
        AN(n)=AN(n)+2*p(m)*cos (2*pi*t(m)*n)*dt/T;
        BN(n)=BN(n)+2*p(m)*sin (2*pi*t(m)*n)*dt/T;
    end
end
end

function [A0,MN,PHIN] = FourierCoefficients2(p,dt,N,T)
[A0,AN,BN] = FourierCoefficients(p,dt,N,T);
MN=zeros(N,1);
PHIN=zeros(N,1);
for n=1:N
    MN(n)=sqrt(AN(n)^2+BN(n)^2);
    PHIN(n)=angle(AN(n)+i*BN(n));
end
end