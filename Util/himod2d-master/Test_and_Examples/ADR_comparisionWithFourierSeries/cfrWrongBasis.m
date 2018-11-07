%%%
% Matteo Aletti, teo.aletti@gmail.com
%
% Comparison with standard Fourier series.
% created 20 novembre 2014
%%%

close all
clear all

addpath ../../Core
addpath ../../Util

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%            Dominio                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nd = 1; %number of sub-domains
L =1; 
domain=[0,L]; %x-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%       Soluzione Esatta                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test=3; % 3 = Educated Basis, 4 = Fourier Series
switch test
%     case {1,2}
%         Ues=sym('x*y+x+y+exp(2*x*y-y)-1 - y*y*(2*x + exp(2*x - 1))');
    case {3,4}
        Ues=sym('x*y+x+y+exp(2*x*y-y)-1 - y*y*(2*x + exp(2*x - 1)+1/10*(2*exp(2*x - 1)*(2*x - 1) - 4*exp(2*x - 1) - 6*x + 2))');
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficienti problema %
%%%%%%%%%%%%%%%%%%%%%%%%%

s            = 2;
bx           = 20;
by           = 0;
mu           = 1;
coeff_robin  = 3;

verb = 1;
expsol = false;

if(verb)
    display(Ues)
    Uf=matlabFunction(Ues);
    x=linspace(0,L);
    y=linspace(0,1);
    [X,Y]=meshgrid(x,y);
    figure;
    contourf(X,Y,Uf(X,Y));
    title(['U_{esatta}(x,y)=',char(Ues)])
    if expsol
        export_py_fun( Uf, domain, 0.01, ['esatta',num2str(test)] );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcolo derivate soluzione esatta%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uesy=diff(Ues,'y');% Derivata in direzione x
Uesx=diff(Ues,'x');% Derivata in direzione y
Ueslapl=diff(Ues,2,'y')+diff(Ues,2,'x');% Laplaciano
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inizializzazione simboli%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=sym('S');
Bx=sym('Bx');
By=sym('By');
Mu=sym('Mu');

% From the exact solution, computation of 
% the force term. First symbolically than it is converted
% in a matlab function to be used as force term.
f=(-Mu*Ueslapl+Bx*Uesx+By*Uesy+S*Ues);
f=subs(subs(subs(subs(f,'S',s),'Bx',bx),'By',by),'Mu',mu);
F=matlabFunction(f);

%Controllo dati al bordo
if(verb)
    Robin_0=coeff_robin*subs(Ues,'y',0)-mu*subs(Uesy,'y',0);
    Robin_1=coeff_robin*subs(Ues,'y',1)+mu*subs(Uesy,'y',1);
    Neumann_1=subs(Uesy,'y',1);
    Neumann_0=subs(Uesy,'y',0);
    Dirichlet_0=subs(Ues,'y',0);
    Dirichlet_1=subs(Ues,'y',1);
    
    display(Dirichlet_0)
    display(Neumann_0)
    display(Robin_0)
    
    display(Dirichlet_1)
    display(Neumann_1)
    display(Robin_1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Is u'' compatible with BC?       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(verb)
    Uesyyy=diff(Ueslapl,'y');% u_yyy
    LRobin_0=coeff_robin*subs(Ueslapl,'y',0)-mu*subs(Uesyyy,'y',0);
    LRobin_1=coeff_robin*subs(Ueslapl,'y',1)+mu*subs(Uesyyy,'y',1);
    LDirichlet_0=subs(Ueslapl,'y',0);
    LDirichlet_1=subs(Ueslapl,'y',1);
    display('Here the value of the boundary conditions computed on u_yy, is it compatible?')
    display(LRobin_0)
    display(LRobin_1)
    display(LDirichlet_0)
    display(LDirichlet_1)
end


Neumann_out=subs(Uesx,'x',L);


% Here the inflow dirichlet data and the 
% outflow neumann data are computed and converted
uin=subs(Ues,'x',0);
Uin=matlabFunction(uin);
Nout=matlabFunction(Neumann_out);

%**********************************************%
% BC laterali:                                 %
%        * 'rob': mu du/dnu + chi u = dato     %
%        * 'dir': u=dato		               %
% Attenzione: funziona bene solo nel caso si   %
% utilizzino le stesse condizioni di bordo     %
% laterali su tutti i domini                   %
%**********************************************%
switch(test)
%     case {1}
%         BC_laterali_up='dir';
%         BC_laterali_down='rob';
%         Dirichlet_1=double(subs(Ues,'y',1));
%         Robin_0=double(coeff_robin*subs(Ues,'y',0)-mu*subs(Uesy,'y',0));
%         bc_up={BC_laterali_up};
%         dato_up={Dirichlet_1};
%         bc_down={BC_laterali_down};
%         dato_down={Robin_0};
%     case {2}
%         BC_laterali_up='rob';
%         BC_laterali_down='rob';
%         Dirichlet_1=double(subs(Ues,'y',1));
%         Robin_0=double(coeff_robin*subs(Ues,'y',0)-mu*subs(Uesy,'y',0));
%         bc_up={BC_laterali_up};
%         dato_up={Dirichlet_1};
%         bc_down={BC_laterali_down};
%         dato_down={Robin_0};
    case {3}
        BC_laterali_up='rob';
        BC_laterali_down='rob';
        Robin_1=0;
        Robin_0=0;
        bc_up={BC_laterali_up};
        dato_up={Robin_1};
        bc_down={BC_laterali_down};
        dato_down={Robin_0};
    case {4}
        BC_laterali_up='per';
        BC_laterali_down='per';
        bc_up={BC_laterali_up};
        bc_down={BC_laterali_down};
        
        % only for compatibility
        data_1=0;
        data_0=0;
        dato_up={data_1};
        dato_down={data_0};
end

%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%
muf    = @(x,y) (  mu + 0*x+0*y ); % Diffusione
beta1f = @(x,y) (  bx + 0*x+0*y ); % Trasporto orizzontale
beta2f = @(x,y) (  by + 0*x+0*y ); % Trasporto verticale
sigmaf = @(x,y) (  s + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',muf,'beta1',beta1f,'beta2',beta2f, ...
    'sigma',sigmaf,'coeffrobin',coeff_robin);

dato_neu=project(Nout,32,bc_up{1},bc_down{1},Coeff_forma);
%**********************************************%
Dati = struct('dato_dir',Uin,'force',F,'dato_neu',dato_neu);
%**********************************************%


%**********************************************%
% Dati relativi alla forma del dominio         %
% ATTENZIONE: il codice funziona solo per L=1  %
% e psi_x=0 in realt?? per L!=1 e psi_x!=0      %
% vanno ricontrollati tutti i coefficienti.    %
% Soprattutto in assembla_x                    %
%**********************************************%
a     = @(x)     0.*x; % Equazione del lato inferiore del canale
L     = @(x) 1 + 0.*x; % Spessore del canale
psi_x = @(x)     0.*x; % Riscalatura
Dati_geometrici = struct('a',a,'L',L,'psi_x',psi_x);
%***********************************************%
%                    SOLVER                     %
%***********************************************%
listOfModes=[1 2 4 8 16];
numberOfDifferenteXStep=1;
numberOfDifferentModes=length(listOfModes); %attenzione, dopo m=22 circa diventa instabile numericamente con 32 nodi di quadratura

% Inizializzazioni
listOfHValues=cell(numberOfDifferenteXStep);
mout=zeros(numberOfDifferentModes,numberOfDifferenteXStep);
L2=zeros(numberOfDifferentModes,numberOfDifferenteXStep);
H1=L2;

hx=0.00125;
for i=1:numberOfDifferenteXStep
    listOfHValues{i}=num2str(hx);
    for k=1:numberOfDifferentModes
        m=listOfModes(k);
        mout(k,i)=m;
        fprintf(1,'h=%f, m=%d, ...   ',hx,m);
        [u,a_ril,b_ril]=solver_DD(domain,m,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,-1,true,'');
        title(['h=',listOfHValues{i},', m=',num2str(m)]);
        %[L2(k,i), H1(k,i)] = compute_error(m,domain(1),domain(2),hx, u{1},a_ril,b_ril,bc_up{1}, bc_down{1}, Coeff_forma,Ues,Uesx,Uesy);
    end
    hx=hx/2;
end

%filename=['caso',num2str(test)];
% save(filename,'h','L2','H1')
% save(filename)