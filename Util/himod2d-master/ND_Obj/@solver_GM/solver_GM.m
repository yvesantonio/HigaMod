%%  Classe che implementa un metodo di tipo Geometrical Multiscale.
%
% 
%
% Metodi Pubblici:
%
%         obj=solver_GM() Costruttore
%         Assemble(obj)
%         CreateMesh(obj)
%         ImportMesh(obj)
%         SetMetodo(obj,str)
%         SetCaseTest(obj,casetest)
%         SetTemp(obj,t0,d,df,theta)
%         SetDD(obj,toll,maxiter,thetaDD)
%         SetFisic(obj,Chi,kappa,Cin,Cest,mu)
%         SetGraphics(obj,treD,frontale,cinf,csup,az,el)
%         Solve(obj)
%         PlotSolution(obj,treD)
%         PlotNiter(obj)
%         PlotErroriDD(obj)
%         SetDefault(obj)
%         Restart(obj)
%         MakeMoreImg(obj,k,varargin)
%         MakeVideo(obj)
%
% Metodi Statici:
%         SetPath
%         Clean
%
% Membri Pubblici:
%         n_iter;
%         Solution;
%
classdef solver_GM<handle
   methods(Access=public)
        function obj=solver_GM()
            %   function obj=solver_GM()
            %   E' il costruttore della classe solver_GM
            %   Setta il numero di passi temporali, inizializza a 0
            %   il numero di iterazioni di DD
            %   Setta Cest=0.1*Cin
            %         u1dold=0;
            %   Lascia le altre variabili al default:
			%         timedipendent=true;
            %         dati fisici
            %
            %         Chi=0.1;    Chi=0 vuol dire neumann omogeneo sui lati
            %         kappa=12;
            %         Cin=1;      Quale unita' di misura??
            %         mu=0.05;    Valore sensato: 7.6e-4;
            %
            %         dati geometrici
            %
            %         lun2d=5;
            %         lun1d=10;
            %         larghezza=0.5 mm
            %         spes=0.01;
            %
            %         variabili relative ai diversi metodi
            %
            %         metodo='full2d';
            %         full2d=1;
            %         krpo=0;
            %         casetest=1;
            %
            %         dati loop spaziale
            %
            %         toll=1e-4;
            %         maxiter=100;
            %         thetaDD=0.5;
            %
            %         variabili relative al loop temporale
            %
            %         t_0=0;
            %         dt=0.5;
            %         t_f=25;
            %         t=0;
            %         theta=2./3;
            %
            %         utility per il problema 1D
            %
            %         N=15;
            %         r = 1
            %         nodi;         Nodi della griglia 1d
            %         griglia;
            %         base;
            %
            %         Matrici per il problema 1D
            %
            %         M
            %         C
            %         K
            %         S
            %         bv
            %
            %         Dati griglie, per il post processing
            %
            %         treD=1;
            %         frontale=0;
            %         az=30;                    Angolo visuale orizzontale
            %         el=60;                      Angolo visuale verticale
            %
            %         variabili ausiliarie per il loop
            %
            %         u1dold;
            %         val_neum_old=0;
            %         first=1;
            %         val_neum_interfaccia;
            %         u1dkmeno1;
            %
            %         variabili per settare il tipo di controllo per la convergenza di
            %         DD
            %
            %         dom1d=0;
            %         dom2d=0;
            %         front1d=1;
            %         front2d=1;
            %
            %         variabili di stato
            %
            %         solved=false;
            %         assembled=false;
            %
            %           raccoltaerrorri per ciclo DD
            %       
            %         raccoltaerr;
            %         raccoltaerr1;
            %         raccoltaerr2;
            %         raccoltaerr3;
            %         raccoltaerr4;
            obj.n_passi=ceil((obj.t_f-obj.t_0)/obj.dt)+1;
            obj.n_iter=zeros(obj.n_passi,1);
            obj.Cest=0.1*obj.Cin;
            obj.u1dold=zeros(obj.N+1,1);
            obj.SetGlobal;
            obj.levels=linspace(obj.cinf,obj.csup,10);
        end
        Assemble(obj)
        CreateMesh(obj)
        ImportMesh(obj)
        SetMetodo(obj,str)
        SetCaseTest(obj,casetest)
        SetTemp(obj,t0,d,df,theta,mtimedipendent)
        SetDD(obj,toll,maxiter,thetaDD)
        SetFisic(obj,Chi,kappa,Cin,Cest,mu)
        SetGraphics(obj,treD,frontale,az,el,varargin)
        Solve(obj)
        PlotSolution(obj,treD)
        PlotNiter(obj)
        PlotErroriDD(obj)
        MakeVideo(obj)
        MakeMoreImg(obj,k,varargin)
        SetDefault(obj)
        Restart(obj)
    end
    
    methods(Access=private,Hidden=true)
        %Exporting
        ExportGeom(obj)
        ExportFirst(obj)
        ExportNeum(obj)
        ExportFisic(obj)
        ExportAll(obj)
        Import_profile(obj,num)
        %To update the value of global variable
        SetGlobal(obj)
        %To save the solution
        SaveTimeStep(obj,num)
        %To plot the solution
        PlotSingleTimeStep(obj,num)
        PlotErroriSingleStimeStep(obj,num)
        MakeImg(obj,k,varargin)
		Solver_stazionario(obj)
        Plot_1d_krpo(obj,num)
    end
    
    methods(Static)
        SetPath
        Clean
    end
    
    properties(SetAccess=public)
        %Numero iterazioni di Domain Decomposition
        n_iter;
        %dati post processing
        Solution;
    end
    properties(SetAccess=public)
		%variabile che indica se è non-stazionario o no.
		timedipendent=true;
        %Coefficiente condizione di Robin
        Chi=3;%0.25; chi=0 vuol dire neumann omogeneo sui lati
		%Costante ncessaria alla creazione di alcuni profili di default, si veda (Vele-Villa)
        kappa=12;
		%Concentrazione in ingresso
        Cin=1; %ma la misura??
		%Concentrazione esterna al vaso
        Cest=0.05;
		%Coefficiente di diffusione dell'ossigeno nel sangue
        mu=0.05;%7.6e-4;
        %Lunghezza della parte di dominio coperta dal modello 2D
        lun2d=15;                  % attenzione deve coincidere nel file ND_mesh.edp
		%Lunghezza della parte di dominio coperta dal modello 1D
        lun1d=6;                 % attenzione deve coincidere nel file ND_mesh.edp
		%Larghezza del vaso
        larghezza=1;            %in millimetri
		%Spessore della zona 2D (a scopo solo grafico)
        spes=0.01;
        %Il metodo scelto per risolvere full2d* *krpo* *naive*
        metodo='full2d';
        full2d=1;
        krpo=0;
        casetest=1;
        %dati loop spaziale
        toll=1e-4;
        maxiter=100;
        thetaDD=0.5;
        %variabili relative al loop temporale
        t_0=0;
        dt=0.5;
        t_f=25;
        t=0;
        theta=2./3;
        n_passi;                                    %Numero passi temporali
        %utility per il problema 1D
        N=300;
        r = 1
        nodi;                                     %Nodi della griglia 1d
        griglia;
        base;
        %Matrici per il problema 1D
        M
        C
        K
        S
        bv
        %Dati griglie, per il post processing
        levels
        points2d
        seg2d
        tri2d
        pointsfull2d
        segfull2d
        trifull2d
        points1d
        seg1d
        tri1d
        treD=1;
        frontale=0;
        cinf=0;     %limite inferiore per la scala in Z
        csup=0.12;   %Limite superiore per la scala in Z
        az=30;   %Angolo visuale orizzontale
        el=60;     %Angolo visuale verticale
        %variabile ausiliaria per il loop
        u1dold;
        val_neum_old=0;
        first=1;
        val_neum_interfaccia;
        u1dkmeno1;
        %variabili per settare il tipo di controllo per la convergenza di
        %DD
        dom1d=0;
        dom2d=0;
        front1d=1;
        front2d=1;
        %ausiliaria per il plot
        u1d;
        currentfig;
        %variabili per il controllo
        solved=false;
        assembled=false;
        %controllo errori DD
        raccoltaerr;
        raccoltaerr1;
        raccoltaerr2;
        raccoltaerr3;
        raccoltaerr4;
        %profilo interfaccia
        profile;        
    end
end
