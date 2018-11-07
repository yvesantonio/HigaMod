close all
clear all
clc
addpath ./../../LIBRERIA_UTENTE/
addpath ./../../UTILITA/

%grado elemento finito
grado=2;
%numero elementi
N=20;
%nodi!=vertici degli elementi ~gdl
nodi = [0 : 1/(grado*N) : 1];
%creazione griglia e base
griglia=struttura_griglia(nodi,grado);
base=struttura_base(grado);

K=termine_diffusione(griglia,base,'co_mu');
bv = termine_noto(griglia,base,'co_f');

dati_bordo = struct( ...
  'bc'    , [2 1] , ... % tipo delle condizioni
  'gamma' , [0 0] , ... % parametro gamma
  'r'     , [0 0] ...  % dato
);

[R r] = termine_bordo(griglia,base,dati_bordo);

A=K+R;
b=bv+r;
u_inc=1:1:griglia.dim-1;
u_note=griglia.dim;


A11 = A(u_inc,u_inc);
A12 = A(u_inc,u_note     );
A21 = A(u_note     ,u_inc);
A22 = A(u_note     ,u_note     );

b1 = b(u_inc);
x2 = [-3];  % valori noti di u al bordo

x1 = A11\(b1-A12*x2);

x(u_inc,1) = x1;
x(u_note     ,1) = x2;

[x_plot y_plot] = visualizza_soluzione(griglia,base,x,20);
