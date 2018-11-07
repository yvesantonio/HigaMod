function u=exactsol(y,t)
global nu w
kmax=100;
u=zeros(size(y));
for k=0:kmax
    numer=w*exp(-nu*(2*k+1)^2*pi^2*t)-w*cos(w*t)+nu*(2*k+1)^2*pi^2*sin(w*t);
    denom=(nu^2*(2*k+1)^4*pi^4+w^2)*(2*k+1)*pi;
    gamma=  8 * nu * numer / denom;   
    u=u+gamma*sin((2*k+1)*pi*y);
end