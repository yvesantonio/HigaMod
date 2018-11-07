close all
clear variables

addpath ../../Core
addpath ../../Util



hx=0.01;
%%%
mL=8;
a     = @(x)       0.*x; % Equazione del lato inferiore del canale
L     = @(x) 1. + 0.*x; % Spessore del canale
Dati_geometriciL = struct('a',a,'L',L);
Omega1DL=[0,1]; % omega 1d before step
%%%
mR=20;
a     = @(x) -1.  + 0.*x; % Equazione del lato inferiore del canale
L     = @(x) 2.   + 0.*x; % Spessore del canale
Dati_geometriciR = struct('a',a,'L',L);
Omega1DR=[1,2]; % omega 2d after  step
%%%



mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) (  20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione
chi = 1;
Coeff_forma = struct('mu',mu,'beta1',beta1,'beta2',beta2, ...
    'sigma',sigma,'coeffrobin',chi);
force = @(x,y) 10*(  ((x-1.5).^2+.4*(y-0.25).^2<0.01) + ((x-1.5).^2+.4*(y-0.75).^2<0.01) );
dato_dir=@(y) y.*(1-y);
dato_neu=zeros(mL,1);
DatiL = struct('dato_dir',dato_dir,'force',force,'dato_neu',dato_neu);
dato_dir=zeros(mR,1);
DatiR = struct('dato_dir',dato_dir,'force',force);


dataUp={0};
dataDownL={0};
dataDownR={0};
bc_up={'dir'};
bc_downL={'dir'};
bc_downR={'dir'};


% iteration stuff
residual = 1;
toll = 1.e-5;
MBR=cell(1,1);

c=1;
while( residual > toll )
    display(c)
    %Solve on the left
    [uL,aRilL,bRilL]=solver_DD(Omega1DL,mL,hx, bc_up, bc_downL,dataUp,dataDownL,DatiL,Coeff_forma,Dati_geometriciL,-1,true,'');
    %compute dirichlet data for the right
    dirOld=DatiR.dato_dir;
    DatiR.dato_dir=BFS_interface(uL,aRilL,bRilL,mL,mR,bc_up{1},bc_downL{1},bc_downR{1},Coeff_forma);
    %right part
    [uR,aRilR,bRilR]=solver_DD(Omega1DR,mR,hx, bc_up, bc_downR,dataUp,dataDownR,...
        DatiR,Coeff_forma,Dati_geometriciR,-1,true,'');
    %compute the neumann data for the left
    neuOld=DatiL.dato_neu;
    DatiL.dato_neu=BFS_interfaceNeu(uR,mL,mR,bc_up{1},bc_downL{1},bc_downR{1},Coeff_forma,hx);
    %residual computation
    residual = norm(neuOld-DatiL.dato_neu,'inf')/norm(neuOld,'inf') + ...
        norm(dirOld-DatiR.dato_dir,'inf')/norm(dirOld,'inf');
    %plot
    plot_solution_DD_BFS(mR,aRilR,bRilR,Omega1DR,hx,uR{1},bc_up{1},bc_downR{1},Coeff_forma,...
        mL,aRilL,bRilL,Omega1DL,uL{1},bc_downL{1});
    pause(0.001);
    c=c+1;
    display(residual)
end