function  SetDefault(obj)
%   function  SetDefault(obj)
%
%   Risetta a default moltissime variabili
%
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem
%         %dati fisici
%         obj.Chi=0.1;%0.25; chi=0 vuol dire neumann omogeneo sui lati
%         obj.kappa=12;
%         obj.Cin=1; %ma la misura??
%         obj.Cest=0.1*obj.Cin;
%         obj.mu=0.05;%7.6e-4;
%         %dati geometrici
%         obj.lun2d=5;                  
%         obj.lun1d=10;                 
%         obj.larghezza=0.5;            %in millimetri
%         obj.spes=0.01;
%         %dati loop temporale
%         obj.n_iter=NaN;
%         obj.t_0=0;
%         obj.dt=0.5;
%         obj.t_f=25;
%         obj.theta=2./3;
%         %dati loop spaziale
%         obj.toll=1e-4;
%         obj.maxiter=100;
%         obj.metodo='full2d';
%         obj.thetaDD=0.5;
%         %dati problema 1D
%         obj.N=15;   
%         obj.casetest=1;
%         %dati iniziali
%         obj.u1dold=NaN;
%         obj.val_neum_old=0;        
%         %dati post processing
%         obj.treD=1;
%         obj.Solution=NaN;
%         obj.frontale=0;
%         obj.cinf=0;     %limite inferiore per la scala in Z
%         obj.csup=1.2;   %Limite superiore per la scala in Z
%         obj.az=30;   %Angolo visuale orizzontale
%         obj.el=60;     %Angolo visuale verticale
%         
%         %variabili relative ai diversi metodi
%         obj.full2d=1;
%         obj.krpo=0;
%         obj.w=NaN;
%         %variabili relative al loop temporale
%         obj.t=0;
%         obj.n_passi=NaN;                                    %Numero passi temporali
%         %utility per il problema 1D
%         obj.r = 1;                                                       
%         obj.nodi=NaN;                                       %Nodi della griglia 1d
%         obj.griglia=NaN;
%         obj.base=NaN;
%         %Matrici per il problema 1D
%         obj.M=NaN;
%         obj.C=NaN;
%         obj.K=NaN;
%         obj.S=NaN;
%         obj.bv=NaN;
%         %Dati griglie, per il post processing
%         obj.points2d=NaN;
%         obj.seg2d=NaN;
%         obj.tri2d=NaN;
%         obj.pointsfull2d=NaN;
%         obj.segfull2d=NaN;
%         obj.trifull2d=NaN;
%         obj.points1d=NaN;
%         obj.seg1d=NaN;
%         obj.tri1d=NaN;
%         %variabile ausiliaria per il loop
%         obj.first=1;
%         obj.val_neum_interfaccia=NaN;
%         obj.u1dkmeno1=NaN;
%         %variabili per settare il tipo di controllo per la convergenza di
%         %DD
%         obj.dom1d=0;
%         obj.dom2d=0;
%         obj.front1d=1;
%         obj.front2d=1;
%         %ausiliaria per il plot
%         obj.u1d=NaN;
        %dati fisici
		obj.timedipendent=true;
        obj.Chi=0.1;%0.25; chi=0 vuol dire neumann omogeneo sui lati
        obj.kappa=12;
        obj.Cin=1; %ma la misura??
        obj.Cest=0.1*obj.Cin;
        obj.mu=0.05;%7.6e-4;
        %dati geometrici
        obj.lun2d=5;                  
        obj.lun1d=10;                 
        obj.larghezza=0.5;            %in millimetri
        obj.spes=0.01;
        %dati loop temporale
        obj.n_iter=NaN;
        obj.t_0=0;
        obj.dt=0.5;
        obj.t_f=25;
        obj.theta=2./3;
        %dati loop spaziale
        obj.toll=1e-4;
        obj.maxiter=100;
        obj.metodo='full2d';
        obj.thetaDD=0.5;
        %dati problema 1D
        obj.N=15;   
        obj.casetest=1;
        %dati iniziali
        obj.u1dold=NaN;
        obj.val_neum_old=0;        
        %dati post processing
        obj.treD=1;
        obj.Solution=NaN;
        obj.frontale=0;
        obj.cinf=0;     %limite inferiore per la scala in Z
        obj.csup=1.2;   %Limite superiore per la scala in Z
        obj.az=30;   %Angolo visuale orizzontale
        obj.el=60;     %Angolo visuale verticale
        
        %variabili relative ai diversi metodi
        obj.full2d=1;
        obj.krpo=0;
        obj.w=NaN;
        %variabili relative al loop temporale
        obj.t=0;
        obj.n_passi=NaN;                                    %Numero passi temporali
        %utility per il problema 1D
        obj.r = 1;                                                       
        obj.nodi=NaN;                                       %Nodi della griglia 1d
        obj.griglia=NaN;
        obj.base=NaN;
        %Matrici per il problema 1D
        obj.M=NaN;
        obj.C=NaN;
        obj.K=NaN;
        obj.S=NaN;
        obj.bv=NaN;
        %Dati griglie, per il post processing
        obj.points2d=NaN;
        obj.seg2d=NaN;
        obj.tri2d=NaN;
        obj.pointsfull2d=NaN;
        obj.segfull2d=NaN;
        obj.trifull2d=NaN;
        obj.points1d=NaN;
        obj.seg1d=NaN;
        obj.tri1d=NaN;
        %variabile ausiliaria per il loop
        obj.first=1;
        obj.val_neum_interfaccia=NaN;
        obj.u1dkmeno1=NaN;
        %variabili per settare il tipo di controllo per la convergenza di
        %DD
        obj.dom1d=0;
        obj.dom2d=0;
        obj.front1d=1;
        obj.front2d=1;
        %ausiliaria per il plot
        obj.u1d=NaN;
        obj.ExportAll;
        obj.SetGlobal;
        obj.assembled=false;
end
