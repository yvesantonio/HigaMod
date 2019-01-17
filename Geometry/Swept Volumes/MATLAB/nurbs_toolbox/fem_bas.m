function [femb,femb_x] = fem_bas(n,x,xl,xr)

  xr = 1;
  xl = 0;
  h = (xr-xl);
  xm = (xl+xr)/2;

  if (n > 2) 

    disp('Error: degree of the fem basis is too large');

  elseif (n < 1) 

    disp('Error: degree of the fem basis is too small');

  % Linear Finite Element

  elseif (n == 1)

    femb = [(xr-x)/h; (x-xl)/h];
    femb_x = [-1/h;1/h];

  % Quadratic Finite Element

  elseif (n == 2)

    femb = [-(x-xr)*(x-xm)/(h*(xl-xm)); (x-xl)*(x-xr)/((xm-xl)*(xm-xr)); (x-xl)*(x-xm)/(h*(xr-xm))];  
    femb_x = [-((2*x-xr-xm))/(h*(xl-xm)); ((2*x-xl-xr))/((xm-xl)*(xm-xr)); ((2*x-xl-xm))/(h*(xr-xm))];

  end; 
  return
end