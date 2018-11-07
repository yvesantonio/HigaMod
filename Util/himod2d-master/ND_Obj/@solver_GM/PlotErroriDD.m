function PlotErroriDD(obj)
% function PlotErroriDD(obj)
% Lancia per ogni istante di tempo la funzione PlotErroriSingleTimeStep
for i=1:size(obj.Solution,1)
    obj.PlotErroriSingleTimeStep(i)
end
end