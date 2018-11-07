function boundary=interp_inflow(ril,mb,boundary,pesi,dato_dir)
%
%   function boundary=interp_inflow(ril,mb,boundary,pesi)
%
%   ril gi?? nei nodi di quadratura
%   TODO a cosa serve passare boundary?
%  
%
m=size(mb,2);
for i = 1 : m
        % Tolgo il rilevamento u_0 = u - l all' INFLOW
        boundary(i) = integrate( ((dato_dir-ril).*mb(:,i))', pesi ); %ATTENZIONE manca jacobiano (L)
end