function []=SetCaseTest(obj,casetest)
%   function []=SetCaseTest(obj,casetest)
%
%   Cambia i dati iniziali e, in alcuni casi, alcune variabili fisiche o
%   geometriche. Per ora sono disponibili valori da 1 a 3, 1 e 2 sono stati
%   controllati.
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem

obj.casetest=casetest;
switch casetest
    case 1
        obj.u1dold=obj.Cest*ones(obj.N+1,1);%Dato di Cauchy (iniziale)
        obj.val_neum_old=0;
    case 2
        obj.u1dold=obj.Cest*ones(obj.N+1,1);%0.9*obj.Cin*ones(obj.N+1,1);%da correggere con quad(vele-villa) invece che 0.9
        obj.val_neum_old=0;
    case 3
        obj.u1dold=obj.Cest*ones(obj.N+1,1);
        obj.val_neum_old=0;
    otherwise
        display('In SetCaseTest: errore switch casetest non riconosciuto, matlab');
end
obj.SetGlobal;
obj.ExportAll;
obj.assembled=false;
end
