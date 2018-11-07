%%
%
% Matteo Aletti, teo.aletti@gmail.com
%
%
% Verifica ordine di convergenza.
% created Ottobre 2014
%
%%

close all
clear most

addpath ../../../Core
addpath ../../../Util
addpath ../../../MapHandler

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

%
% Ues=sym('3*y*(1-y) + x*x*4*y^2*(1-y)*(0.75+8*x^2*y+8*x*y^2)*(x-1)^2+(1-y)^2');

% f=f(x,y); t.c. f_y|y=1 = -1 
% f=x^2/2 -x^3/3 - y
% f_y= -1
% f_x=x*(1-x)
% Ues=sym('exp(f)');
% rob=exp(f)(1+f_y)
% neux=exp(f)(1+f_x)

%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficienti problema %
%%%%%%%%%%%%%%%%%%%%%%%%%

s=2;
bx=20;
by=3;
mu=1;
coeff_robin=1;

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
        export_py_fun(Uf,cutx,0.01,['esatta_legendre']);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calcolo derivate soluzione esatta%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uesy=diff(Ues,'y')% Derivata in direzione x
Uesx=diff(Ues,'x')% Derivata in direzione y
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
% if(verb)
%     Uesyyy=diff(Ueslapl,'y');% u_yyy
%     LRobin_0=coeff_robin*subs(Ueslapl,'y',0)-mu*subs(Uesyyy,'y',0);
%     LRobin_1=coeff_robin*subs(Ueslapl,'y',1)+mu*subs(Uesyyy,'y',1);
%     LDirichlet_0=subs(Ueslapl,'y',0);
%     LDirichlet_1=subs(Ueslapl,'y',1);
%     display('Here the value of the boundary conditions computed on u_yy, is it compatible?')
%     display(LRobin_0)
%     display(LRobin_1)
%     display(LDirichlet_0)
%     display(LDirichlet_1)
% end

if(verb)
    Neumann_out=subs(Uesx,'x',L);
    display('Value of u_x on the outflow, should be equal to zero: ')
    display(Neumann_out);
    Neumann_in=subs(Uesx,'x',0);
    display('Value of u_x on the inflow, should be equal to zero: ')
    display(Neumann_in);
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

BC_laterali_up='rob';
BC_laterali_down='dir';
Dirichlet_0=double(subs(Ues,'y',0));
Robin_1=double(coeff_robin*subs(Ues,'y',1)+mu*subs(Uesy,'y',1));
bc_up={BC_laterali_up};
dato_up={Robin_1};
bc_down={BC_laterali_down};
dato_down={Dirichlet_0};

%**********************************************%
Dati = struct('dato_dir',Uin,'force',F);
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

% %***********************************************%
% %                    SOLVER                     %
% %***********************************************%
% m_=[1 2 4 8 16 32];
% NH=5;
% NM=length(m_); %attenzione, dopo m=22 circa diventa instabile numericamente con 32 nodi di quadratura
% 
% % Inizializzazioni
% h=cell(NH);
% mout=zeros(NM,NH);
% L2=zeros(NM,NH);
% H1=L2;
% haty=sym('y');
% [Map]=provideMap(haty);
% hx=0.2;
% for i=1:NH
%     hx=hx/2;
%     h{i}=num2str(hx);
%     mesh_x=cutx(1):hx:cutx(2);
%     for k=1:NM
%         m=m_(k);
%         mout(k,i)=m;
%         fprintf(1,'h=%f, m=%d, ...   ',hx,m);
% 
% %       [u]=solver_leg(mesh_x, Coeff_forma, m, F,Map,bc_laterali_up,bc_laterali_down);
%         [u,a_ril,b_ril]=solver_DD(cutx,m,hx, bc_up, bc_down,dato_up,dato_down,Dati,Coeff_forma,Dati_geometrici,-1,false,'');
%         title(['h=',h{i},', m=',num2str(m)]);
%         [L2(k,i), H1(k,i)] = compute_error(m,cutx(1),cutx(2),hx, u{1},a_ril,b_ril,bc_up{1}, bc_down{1}, Coeff_forma,Ues,Uesx,Uesy);
%         
%     end
% end
% 
% filename=['caso',num2str(test)];
% save(filename,'h','L2','H1')
% 
% % Exporter per tikz
% for k=1:NM
%     L2(k,NH+1)=L2(1,1)/(2^k);
%     L2(k,NH+2)=L2(1,1)/(4^k);
%     L2(k,NH+3)=L2(1,1)/(8^k);
%     L2(k,NH+4)=L2(1,1)/(16^k);
%     H1(k,NH+1)=H1(1,1)/(2^k);
%     H1(k,NH+2)=H1(1,1)/(4^k);
%     H1(k,NH+3)=H1(1,1)/(8^k);
%     H1(k,NH+4)=H1(1,1)/(16^k);
% end
% fil=fopen(['Caso',num2str(test),'_L2.dat'],'w');
% fprintf(fil,'#h= ');
% for i=1:NH-1
%     fprintf(fil,[h{i},',']);
% end
% fprintf(fil,[h{NH},'\n']);
% for k=1:NM
%     fprintf(fil,'%d, ',mout(k,1));
%     for i=1:NH+3
%         fprintf(fil,'%12.16f, ',L2(k,i));
%     end
%     fprintf(fil,'%12.16f\n',L2(k,NH+4));
% end
% fprintf(fil,'#-----\n');
% fclose(fil);
% 
% fil=fopen(['Caso',num2str(test),'_H1.dat'],'w');
% fprintf(fil,'#h= ');
% for i=1:NH-1
%     fprintf(fil,[h{i},',']);
% end
% fprintf(fil,[h{NH},'\n']);
% for k=1:NM
%     fprintf(fil,'%d, ',mout(k,1));
%     for i=1:NH+3
%         fprintf(fil,'%12.16f, ',H1(k,i));
%     end
%     fprintf(fil,'%12.16f\n',H1(k,NH+4));
% end
% fprintf(fil,'#-----\n');
% fclose(fil);
