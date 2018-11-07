close all
clear all
addpath ..
global a L psi_x coeffrobin mu beta1 beta2 sigma force

numero_modi = 20;                                    % Numero modi
Estremi     = [0,10];								 % Estremi del dominio
degree_x    = 1;                                     % Grado degli elementi finiti
hx          = 0.1;                                   % Passo spaziali								
chi         = 0.25;

bc_up     = 'dir';
dato_up   = 3;
bc_down   = 'rob';
dato_down = 3;

%forzrobin  = -chi*cest;    % mu*c_nu+coeffrobin*c=forzrobin
coeffrobin = chi;

% TODO Aggiungere proiezione dato.
% Identifica le condizioni al bordo fisico (INFLOW)
BC_l = 'dir';
% Valore della condizione al bordo fisico (INFLOW)
boundary_l = zeros(numero_modi, 1);
% Identifica le condizioni al bordo fisico (OUTFLOW)
BC_r = 'neu';
% Valore della condizione al bordo fisico (OUTFLOW)
boundary_r = zeros(numero_modi,1); % attenzione moltiplicare per mu se serve!!!

% Asse del centro del canale, potrebbe essere una funzione di x
a     = @(x) 0.*x;
% Spessore del canale, potrebbe essere una funzione di x
L     = @(x) 1+ 0.*x;
% Riscalatura
psi_x = @(x)  0.*x;

%Coefficienti della forma bilineare
mu    = @(x,y) (  1 + 0*x+0*y ); % Diffusione
beta1 = @(x,y) ( 20 + 0*x+0*y ); % Trasporto orizzontale
beta2 = @(x,y) (  0 + 0*x+0*y ); % Trasporto verticale
sigma = @(x,y) (  0 + 0*x+0*y ); % Reazione
force = @(x,y) 0;%( ((x-1.5).^2 + (y-0.25).^2) < 0.01) + ( ((x-1.5).^2 + (y-0.75).^2) < 0.01); %forzante

% x=linspace(Estremi(1),Estremi(2),100);
% y=linspace(0,1,100);
% [X,Y]=meshgrid(x,y);
% surf(X,Y,force(X,Y))
% axis equal
% pause

%Solver
[u]=solver(Estremi, numero_modi, degree_x, hx, BC_l, boundary_l, BC_r, boundary_r, bc_up, bc_down, dato_up, dato_down);
axis equal

pause

! FreeFem++ -ne calcola_esatta.edp
