function PlotNiter(obj)
% function PlotNiter(obj)
% Plotta il numero di iterazioni di DD
%

t=obj.dt:obj.dt:(obj.n_passi-1)*obj.dt;
figure(3)

plot(t,obj.n_iter(2:end))
xlabel('t');
ylabel('# iter');
title('numero di iterazioni di Domain-Decomposition')
end