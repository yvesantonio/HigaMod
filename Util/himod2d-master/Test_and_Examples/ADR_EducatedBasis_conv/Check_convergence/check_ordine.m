%%
%
% Matteo Aletti, teo.aletti@gmail.com
%
%
% Verifica ordine di convergenza.
% created 28 aprile 2013
%
%%

close all
clear all

addpath ../../../Core
addpath ../../../Util

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Attenzione, il numero di nodi di quadratura (normalmente 32)
%   pu?? influire sulla convergenza al crescere del numero di modi.
%   Se si vuole indagare il caso con 32 modi occorre fissare almeno 64
%   Per farlo modificare Core/build_system.m
%   Caso test 1:
%              condizioni di Robin non omogenee in alto
%              condizioni di Dirichlet omogenee in basso
%              EH1: leggera superconvergenza (mezzo ordine?)
%			   EL2: leggera superconvergenza (mezzo ordine?)
%              Coefficienti di Fourier da guardare ancora.
%              Attenzione l'ordine di infinitesimo non pu?? essere 2 perch??
%              la funzione non ?? compatibile con i dati al bordo, senza
%              rilevamento.
%   Caso test 2:
%              condizioni di Robin     omogenee in alto
%              condizioni di Dirichlet omogenee in basso
%              EH1: leggera superconvergenza (mezzo ordine?)
%			   EL2: leggera superconvergenza (mezzo ordine?)
%              Coefficienti di Fourier da guardare ancora.
%   Caso test 3:
%		       condizioni di Robin     omogenee in alto
%              condizioni di Dirichlet omogenee in basso
%              EH1: sembra ragionevole, ma non chiaro.
%			   EL2: leggera superconvergenza (mezzo ordine?)
%              Coefficienti di Fourier da guardare ancora. Caso complicato
%              sarebbe il caso di lisciare le m per vedere un po' meglio
%			   Lisciando si osserva che c'?? un chiaro trend, ma ?? necessario infittire
%              la griglia.
%   Caso test 4:
%		       condizioni di Robin     omogenee in alto
%              condizioni di Dirichlet omogenee in basso
%              EH1: lungo plateau poi rispecchia
%			   EL2: lungo plateau poi rispecchia
%              La convergenza cmq si verifica solo dopo 8 modi.
%              Infatti, dai coeff Fourier si osserva che il comportamento
%              asintotico si realizza a partire da circa 10 modi.
%   Caso test 5:
%			   condizioni di Robin     omogenee in alto
%              condizioni di Dirichlet omogenee in basso
%              extracompatibile
%              si osserva extraconvergenza
%              ancora da misurare.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%            Dominio                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nd = 1;

%y \in (0,1)
% Lunghezza
L =1;
% vettore "cutx"
cutx=[0,L];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%       Soluzione Esatta                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expsol=true;
verb=true;
test=10;
switch(test)
    case 1  % Caso non omogeneo caso del paper!
        Ues=sym('4*y^2*(1-y)*(0.75+8*x^2*y+8*x*y^2)*(x-1)^2+(1-y)^2');
    case 2  % Omogeneo
        Ues=sym('y^2*(1-y)^2*exp(x)*(x-1)^2');
    case 3
        Ues=sym('y^2*(1-y)^3*(y-1/2)*exp(10*x)*(x-1)^2');
    case 4
        Ues=sym('(1-y)*(y-1/2)*exp(10*x*y^2)*y^2*(x-1)^2');
    case 5  % caso super-convergence (voluto) paper case
        Ues=sym('y^4*(1-y)^4*exp(x)*(x-1)^2');
    case 6  % Robin Robin super convergence (why??)
        Ues=sym('3*y*(1-y)+3+4*y^2*(1-y)^2*(0.75+8*x^2*y+8*x*y^2)*(x-1)^2+(1-y)^2+y-2');
    case 7  % Robin Robin super convergence (why???) attenzione una frequenza ogni due ha coefficiente di fourier nullo
        Ues=sym('3*y*(1-y)*exp(y^2*(1-y)^2*3*(sin(2*pi*x))^2)+3');
    case 8  % Robin Robin super convergece ????
        Ues=sym('y*(1-y)*exp(sin( 20*y^3*(1-y)^2*(x-1)^2 )) + 2');
    case 9  % dir dir
        Ues=sym('y*(1-y)*exp(sin( 20*y^3*(1-y)^2*(x-1)^2 ))');
    case 10  % neu neu
        Ues=sym('y^2*(1-y)^2*exp(sin( 20*y^3*(1-y)^2*(x-1)^2))+1');
    case 11 %rob dir
        Ues=sym('y^2*(1-y)^3*(y-1/2)*exp(10*x)*(x-1)^2');
    case 12
        Ues=sym('10^5*exp( - 12*0.1*(y* - 1/2 ) ^ 2 + (-2.0/3.0/( 1.0 + x)*0.1^3*y^3+1.0/( 1.0 + x)*y^2*0.1^2)*70  )*(1-x)^2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficienti problema %
%%%%%%%%%%%%%%%%%%%%%%%%%
switch(test)
    case {1,2,5,8}
        s=2;
        bx=20;
        by=0;
        mu=1;
        coeff_robin=10000;
    case {3,4}
        s=2;
        bx=20;
        by=0;
        mu=1;
        coeff_robin=3;
    case {6,7,9,10,11}
        s=2;
        bx=20;
        by=2;
        mu=1;
        coeff_robin=1;
    case {12}
        s=1;
        bx=10;
        by=5;
        mu=1;
        coeff_robin=2;
        gamma12=coeff_robin/mu;
        Ues=subs(Ues,'gamma12',gamma12);
end
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
        export_py_fun(Uf,cutx,0.01,['esatta',num2str(test)]);
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

% Calcolo forzante
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
    display(Robin_0)
    display(Robin_1)
    display(Dirichlet_0)
    display(Dirichlet_1)
    display(Neumann_0)
    display(Neumann_1)
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

if(verb)
    Neumann_out=subs(Uesx,'x',L);
    display('Value of u_x on the outflow, should be equal to zero: ')
    display(Neumann_out);
end

% Conversione dato dirichlet inflow
uin=subs(Ues,'x',0);
Uin=matlabFunction(uin);

% Conversione dato dirichlet outflow eventualmente se si impone dirichlet all'outflow
% uout=subs(Ues,'x',L);
% Uout=matlabFunction(uout);

%**********************************************%
% BC laterali:                                 %
%        * 'rob': mu du/dnu + chi u = dato     %
%        * 'dir': u=dato		               %
% Attenzione: funziona bene solo nel caso si   %
% utilizzino le stesse condizioni di bordo     %
% laterali su tutti i domini                   %
%**********************************************%
switch(test)
    case {1,2,3,4,5}
        BC_laterali_up='dir';
        BC_laterali_down='rob';
        Dirichlet_1=double(subs(Ues,'y',1));
        Robin_0=double(coeff_robin*subs(Ues,'y',0)-mu*subs(Uesy,'y',0));
        bc_up={BC_laterali_up};
        dato_up={Dirichlet_1};
        bc_down={BC_laterali_down};
        dato_down={Robin_0};
    case {6,7,8,12}
        BC_laterali_up='rob';
        BC_laterali_down='rob';
        Robin_1=double(coeff_robin*subs(Ues,'y',1)+mu*subs(Uesy,'y',1));
        Robin_0=double(coeff_robin*subs(Ues,'y',0)-mu*subs(Uesy,'y',0));
        bc_up={BC_laterali_up};
        dato_up={Robin_1};
        bc_down={BC_laterali_down};
        dato_down={Robin_0};
    case 9
        BC_laterali_up='dir';
        BC_laterali_down='dir';
        Dirichlet_0=double(subs(Ues,'y',0));
        Dirichlet_1=double(subs(Ues,'y',1));
        bc_up={BC_laterali_up};
        dato_up={Dirichlet_1};
        bc_down={BC_laterali_down};
        dato_down={Dirichlet_0};
    case 10
        BC_laterali_up='neu';
        BC_laterali_down='neu';
        Neumann_0=double(subs(Uesy,'y',0));
        Neumann_1=double(-subs(Uesy,'y',1));
        bc_up={BC_laterali_up};
        dato_up={Neumann_1};
        bc_down={BC_laterali_down};
        dato_down={Neumann_0};
    case 11
        BC_laterali_up='rob';
        BC_laterali_down='dir';
        Dirichlet_0=double(subs(Ues,'y',0));
        Robin_1=double(coeff_robin*subs(Ues,'y',1)+mu*subs(Uesy,'y',1));
        bc_up={BC_laterali_up};
        dato_up={Robin_1};
        bc_down={BC_laterali_down};
        dato_down={Dirichlet_0};
end



%**********************************************%
Dati = struct('dato_dir',Uin,'force',F);
%**********************************************%

%**********************************************%
% Grafico del dato in ingresso                 %
%**********************************************%
%if(verb)
%    figure;
%    fplot(Uin,[0,1]);
%    title('Dato in ingresso');
%    xlabel('y');
%    ylabel('u_{in}');
%end

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
%**********************************************%


%**********************************************%
% Coefficienti della forma bilineare           %
%**********************************************%
muf    = @(x,y) (  mu + 0*x+0*y ); % Diffusione
beta1f = @(x,y) ( bx + 0*x+0*y ); % Trasporto orizzontale
beta2f = @(x,y) (  by + 0*x+0*y ); % Trasporto verticale
sigmaf = @(x,y) (  s + 0*x+0*y ); % Reazione

Coeff_forma = struct('mu',muf,'beta1',beta1f,'beta2',beta2f, ...
    'sigma',sigmaf,'coeffrobin',coeff_robin);
%**********************************************%
% tentativo di stimare dalla soluzione esatta  %
% i suoi coefficienti di fourier     |         %
% under-development                  V         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(verb)
    figure
    soglia=10;
    set(gca,'FontUnits','points','FontSize',12);
    uk=compute_fourier_coefficient(matlabFunction(Ues),Coeff_forma,bc_up{1},bc_down{1});
    [ordine,C]=compute_asimptotyc_order(uk,soglia);
    loglog(uk,'*');
    k=1:length(uk);
    hold on
    loglog(k,C*k.^(-ordine),'--')
    fprintf(1,'Ordine di infinitesimo stimato dei coefficienti di Fourier:\n');
    fprintf(1,'                  %.2f, ottenuto a partire da k=%d\n',ordine,soglia);
    xlabel('k')
    ylabel('u_k')
    title('Coefficienti di Fourier di Ues, comportamento asintotico')
end

%***********************************************%
%                    SOLVER                     %
%***********************************************%
m_=[1 2 4 8 16 32];
%m_=[1 2];
NH=5;
NM=length(m_); %attenzione, dopo m=22 circa diventa instabile numericamente con 32 nodi di quadratura

% Inizializzazioni
h=cell(NH);
mout=zeros(NM,NH);
L2=zeros(NM,NH);
H1=L2;

hx=0.2;
for i=1:NH
    hx=hx/2;
    h{i}=num2str(hx);
    for k=1:NM
        m=m_(k);
        mout(k,i)=m;
        fprintf(1,'h=%f, m=%d, ...   ',hx,m);
        [u,a_ril,b_ril]=solver_DD(cutx,m,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,-1,false,'');
        title(['h=',h{i},', m=',num2str(m)]);
        [L2(k,i), H1(k,i)] = compute_error(m,cutx(1),cutx(2),hx, u{1},a_ril,b_ril,bc_up{1}, bc_down{1}, Coeff_forma,Ues,Uesx,Uesy);
        
    end
end

filename=['caso',num2str(test)];
save(filename,'h','L2','H1')

% Exporter per tikz
for k=1:NM
    L2(k,NH+1)=L2(1,1)/(2^k);
    L2(k,NH+2)=L2(1,1)/(4^k);
    L2(k,NH+3)=L2(1,1)/(8^k);
    L2(k,NH+4)=L2(1,1)/(16^k);
    H1(k,NH+1)=H1(1,1)/(2^k);
    H1(k,NH+2)=H1(1,1)/(4^k);
    H1(k,NH+3)=H1(1,1)/(8^k);
    H1(k,NH+4)=H1(1,1)/(16^k);
end
fil=fopen(['Caso',num2str(test),'_L2.dat'],'w');
fprintf(fil,'#h= ');
for i=1:NH-1
    fprintf(fil,[h{i},',']);
end
fprintf(fil,[h{NH},'\n']);
for k=1:NM
    fprintf(fil,'%d, ',mout(k,1));
    for i=1:NH+3
        fprintf(fil,'%12.16f, ',L2(k,i));
    end
    fprintf(fil,'%12.16f\n',L2(k,NH+4));
end
fprintf(fil,'#-----\n');
fclose(fil);

fil=fopen(['Caso',num2str(test),'_H1.dat'],'w');
fprintf(fil,'#h= ');
for i=1:NH-1
    fprintf(fil,[h{i},',']);
end
fprintf(fil,[h{NH},'\n']);
for k=1:NM
    fprintf(fil,'%d, ',mout(k,1));
    for i=1:NH+3
        fprintf(fil,'%12.16f, ',H1(k,i));
    end
    fprintf(fil,'%12.16f\n',H1(k,NH+4));
end
fprintf(fil,'#-----\n');
fclose(fil);
