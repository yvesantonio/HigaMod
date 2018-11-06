function [p,isoGeoMap,error] = bspevalKalman(d,sol,k,mesh) 
%  
% Function Name: 
%  
%   bspeval - Evaluate a univariate B-Spline. 
%  
% Calling Sequence: 
%  
%   p = bspeval(d,c,k,u) 
%  
% Parameters: 
%  
%   d	: Degree of the B-Spline. 
%  
%   c	: Control Points, matrix of size (dim,nc). 
%  
%   k	: Knot sequence, row vector of size nk. 
%  
%   u	: Parametric evaluation points, row vector of size nu. 
%  
%   p	: Evaluated points, matrix of size (dim,nu) 
%  
% Description: 
%  
%   Evaluate a univariate B-Spline. This function provides an interface to 
%   a toolbox 'C' routine. 

nu = numel(mesh);

[mc,nc] = size(sol);

isoGeoMap = zeros(nc,nu);

p = zeros(mc,nu);

N = zeros(d+1,1);                                
                                                
for col=1:nu                                    
                                                
    s = findspan(nc-1, d, mesh(col), k);
    N = basisfun(s,mesh(col),d,k);    
    
    tmp1 = s - d + 1;
    
    for row=1:mc                                
        tmp2 = 0;                               
        for i=0:d                               
            tmp2 = tmp2 + N(i+1)*sol(row,tmp1+i);
               
            isoGeoMap(tmp1+i,col) = N(i+1);
        end                                     
        p(row,col) = tmp2;                      
    end                                         
end

p_prime = sol * isoGeoMap;
error = abs(p - p_prime); 