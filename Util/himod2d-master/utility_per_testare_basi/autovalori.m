clear all
close all
%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%
chi=4.456;
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) ( 0 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigmav = @(x,y) (  0 + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigmav,'coeffrobin',chi);

addpath ../Core
addpath ../
m=30;
L=2;
bc_up='rob';
bc_down='rob';
muval=mu(0,0)/L;
sigma=chi;
stlv  = @(lambda) 2*muval*lambda + tan(lambda).*(sigma - muval.^2*lambda.^2/sigma);
dstlv=@(x) 0.*x;

lambda=compute_autovalues(stlv,dstlv,m,bc_up,bc_down,Coeff_forma);
lambda/L
k=1:m;
xv=(k-1)*pi;
plot(k,xv,'red','linewidth',3)
hold on
plot(lambda,'.','MarkerSize',15)


dl=lambda(2:end)-lambda(1:end-1);
for i=1:m
stima_c(i)=lambda(i)-i*pi;
end


nome=['./','numerical'];
file=fopen(nome,'w+');
v=1:m;
stringa=['# m       lambda\n'];
fprintf(file,stringa);
for i=1:m
	stringa=[num2str(v(i)),'    ',num2str(lambda(i)),'\n'];
	fprintf(file,stringa);
end
fclose(file);

nome=['./','teorical'];
file=fopen(nome,'w+');
v=1:m;
stringa=['# m       lambda\n'];
fprintf(file,stringa);
for i=1:m
	stringa=[num2str(v(i)),'    ',num2str(xv(i)),'\n'];
	fprintf(file,stringa);
end
fclose(file);