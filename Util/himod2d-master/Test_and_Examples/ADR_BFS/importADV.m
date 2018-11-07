% create a grid of quadrature point on right domain
addpath ../../Core
addpath ../../Util
addpath ./Util
hx=0.01;
nqnx = 4;
nqny = 64;
ne = round((1.-0.)/hx);
nx = ne+1;
[~, xq, wxq] = gausslegendre( nqnx ); % nqnx non viene cambiato
[~, yq, wyq] = gausslegendre( nqny ); % nqny non viene cambiato
meshx = zeros(nx,1);
meshx(1) = 0;
for i=2:nx
    meshx(i) = meshx(i-1)+hx;
end
mesh_xx = zeros( ne*nqnx, 1); % Nodi
for i=1:ne
    [mesh_xx( (i-1)*nqnx+1 : i*nqnx ), ~ ] = quadrature_rule( meshx(i), meshx(i+1), xq, wxq);
end

[x,YQ] = meshgrid( mesh_xx, yq );
L=@(x) 1 +0*x;
a=@(x) 0 +0*x;
L_computed     = L(x)'; % spessore
a_computed     = a(x)'; % cordinate y inferiore del canale

y = a_computed' + L_computed'.*YQ;
XL=x;
YL=y;
% create a grid of quad point on left domain
nqnx = 4;
nqny = 64;
ne = round((2.-1.)/hx);
nx = ne+1;
[~, xq, wxq] = gausslegendre( nqnx ); % nqnx non viene cambiato
[~, yq, wyq] = gausslegendre( nqny ); % nqny non viene cambiato
meshx = zeros(nx,1);
meshx(1) = 1;
for i=2:nx
    meshx(i) = meshx(i-1)+hx;
end
mesh_xx = zeros( ne*nqnx, 1); % Nodi
for i=1:ne
    [mesh_xx( (i-1)*nqnx+1 : i*nqnx ), ~ ] = quadrature_rule( meshx(i), meshx(i+1), xq, wxq);
end

[x,YQ] = meshgrid( mesh_xx, yq );
L=@(x) 2 +0*x;
a=@(x) -1 +0*x;
L_computed     = L(x)'; % spessore
a_computed     = a(x)'; % cordinate y inferiore del canale

y = a_computed' + L_computed'.*YQ;
XR=x;
YR=y;
% save it in a file
mySave(XR,XL,YR,YL);
% launch a FreeFem script to read and evaluate the advection field there
% and save it to a new file
system('/usr/local/bin/FreeFem++ evaluateField.edp');
% open such file and save it as .mat in ./adv/
importfile('./adv/b1.dat'); %b1
importfile('./adv/b2.dat'); %b2

[B1L,B2L,B1R,B2R]=createMatrices(b1,b2,size(XL,1),size(XL,2));

save('./adv/allfields.mat','B1L','B2L','B1R','B2R')
