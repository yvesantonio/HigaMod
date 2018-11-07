function []=Assemble(obj)
%   function []=Assemble(obj)
%   La funzione assembla il problema
%   In particolare:
%   Aggiorna le variabili globali
%   Setta le variabili di ciclo in tempo
%   costruisce le variabili nodi,griglia,base del problema 1D
%
%   Calcola tutte le matrici costanti in tempo (per il problema 1D)
%   e il termine noto.
%   setta il valore di neumann all'interfaccia e il valore di u1dkmeno1
%   Esporta tutti i dati necessari a FreeFem
%   Inizializza il contenitore della soluzione
%   Segna assebled come vero.
    obj.SetGlobal;
	if(obj.timedipendent)
	    obj.n_passi=ceil((obj.t_f-obj.t_0)/obj.dt)+1; %CASO stat non ha senso
	else
		obj.n_passi=1;	
	end
    obj.n_iter=zeros(obj.n_passi,1);
    %inizializzazione problema 1D
    obj.nodi = obj.lun2d : obj.lun1d/(obj.r*obj.N) : obj.lun1d + obj.lun2d;             %Nodi della griglia
    obj.griglia = struttura_griglia(obj.nodi,obj.r);
    obj.base = struttura_base(obj.r);
    obj.M=termine_reazione(obj.griglia,obj.base,'massa');   %Matrice massa assunta costante
    obj.K=termine_diffusione(obj.griglia,obj.base,'coeff_mu'); %Matrice diffusione costante
    obj.bv = termine_noto(obj.griglia,obj.base,'coeff_f'); %Forzante assunta costante
    %inizializzazioni valori iniziali e guess
    obj.val_neum_interfaccia=obj.val_neum_old;  %CASO stat-->??
    obj.u1dkmeno1=obj.u1dold;  %la prima iterata spaziale coincide con quella temporale %CASO stat-->??
    obj.ExportAll;
    if(~obj.krpo)
        obj.S = termine_reazione(obj.griglia,obj.base,'coeff_s'); %Matrice reazione assunta costante
    end
    obj.Solution=cell(obj.n_passi,2); %	Caso stat--> ok
    obj.CreateMesh;
    obj.ImportMesh;
    obj.assembled=true;
end
