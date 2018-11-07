function SetFisic(obj,Chi,kappa,Cin,Cest,mu)
%	function SetFisic(obj,Chi,kappa,Cin,Cest,mu)
%
%   Setta i parametri fisici 
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem
        obj.Chi=Chi;%0.25; chi=0 vuol dire neumann omogeneo sui lati
        obj.kappa=kappa;
        obj.Cin=Cin; %ma la misura??
        obj.Cest=Cest;
        obj.mu=mu;%7.6e-4;
        obj.ExportAll;
        obj.SetGlobal;
        obj.assembled=false;
end
