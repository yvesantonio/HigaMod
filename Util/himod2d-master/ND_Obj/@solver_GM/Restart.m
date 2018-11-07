function Restart(obj)
%   function Restart(obj)
%   risetta t a t_0, solved a false.
%   azzera Solution 
%   ricalcola e reimporta la mesh
%   Reimposta come dati iniziali i dati del casetest memorizzato
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem
    obj.solved=false;
    obj.first=1;
    obj.t=obj.t_0; % CASO stat non ha senso
    obj.ExportAll;
    obj.Solution=cell(obj.n_passi,2); %CASO stat cell(1,2)
    obj.SetCaseTest(obj.casetest); %%CASO stat ??
    obj.CreateMesh;
    obj.ImportMesh;
    obj.SetGlobal;
    obj.assembled=false;    
end
