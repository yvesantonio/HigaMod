function SetDD(obj,toll,maxiter,thetaDD)
%   function SetDD(obj,toll,maxiter,thetaDD)
%
%   Setta la tolleranza per il test di convergenza di DD e anche il suo
%   numero di iterazioni massimo e il parametro thetaDD di rilassamento.
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem
        obj.toll=toll;
        obj.maxiter=maxiter;
        obj.thetaDD=thetaDD;
        obj.ExportAll;
        obj.SetGlobal;
        obj.assembled=false;
end
