function I = integrami(f,x)
%
% function I = integrami(f,x)
% Integra una funzione definita per punti , f e x devono essere due vettori lunghi uguale.
%
I=0;
for i=1:length(x)-1
	dx=x(i+1)-x(i);
	dy=(f(i+1)+f(i))/2;
	I=I+dx*dy;
end
