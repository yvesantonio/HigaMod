function PlotSolution(obj)
% function PlotSolution(obj)
% Lancia per ogni istante di tempo la funzione PlotSingleTimeStep
obj.currentfig=figure;
for i=1:size(obj.Solution,1)
    obj.PlotSingleTimeStep(i)
end
end