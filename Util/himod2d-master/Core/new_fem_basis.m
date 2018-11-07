function [femb,femb_x]=new_fem_basis(degree,meshx, mesh_xx)

% numero di intervalli
  ne = length(meshx)-1;
% Numero di nodi di quadratura in ciascun intervallo
  nqnx = length(mesh_xx) / ne;
% sono solo due righe, gli zeri non vengono mai memorizzati in ogni riga
% c'è un pezzo della funzione di base
  femb = zeros(length(mesh_xx),degree+1);
  femb_x = zeros(length(mesh_xx),degree+1);

  for ie = 1:ne
    xl = meshx(ie);             % lato sinistro dell'intervallo
    xr = meshx(ie+1);           % lato destro dell'intervallo
    for iqn = 1:nqnx   
      i_x = (ie-1)*nqnx+iqn;
      x = mesh_xx(i_x);%estraggo i nodi di integrazione di questo intervallo
      [femb(i_x,:),femb_x(i_x,:)] = fem_bas(degree,x,xl,xr);
    end
  end
return

function [femb,femb_x]=fem_bas(n,x,xl,xr)

  h = (xr-xl);
  xm = (xl+xr)/2;
  if (n > 2) 
    disp('Error: degree of the fem basis is too large');
  elseif (n < 1) 
    disp('Error: degree of the fem basis is too small');
  elseif (n == 1) %ef lineari
    femb = [(xr-x)/h; (x-xl)/h];
    femb_x = [-1/h;1/h];
  elseif (n == 2) %ef quadratici
    femb = [(x-xr)*(x-xm)/(h*(xl-xm)); (x-xl)*(x-xr)/((xm-xl)*(xm-xr)); (x-xl)*(x-xm)/(h*(xr-xm))];
    femb_x = [((2*x-xr-xm))/(h*(xl-xm)); ((2*x-xl-xr))/((xm-xl)*(xm-xr)); ((2*x-xl-xm))/(h*(xr-xm))];
  end;  
return
