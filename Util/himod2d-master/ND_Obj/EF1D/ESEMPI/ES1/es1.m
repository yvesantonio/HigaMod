% Consideriamo il problema:
%
%   -d2u/dx2 + beta*du/dx = 3*x   per x in (0,10)
%   u = 0                         per x=0
%   -du/dx - 3*u = 15             per x=10
%
% essendo beta un parametro reale. La soluzione analitica di tale
% problema è:
%
%   u = C*(1-exp(beta*x)) + 3*x/(2*beta^2)*(beta*x+2).
%
% con
%
%   C = (15+3/beta^2*(160*beta+31))/((beta+3)*exp(10*beta)-3).
%
% In questo esempio, possiamo considerare i due valori beta=1
% (diffusione dominante) e beta=10 (trasporto dominante).

close all
clear all

addpath ./../../LIBRERIA_UTENTE
addpath ./../../UTILITA

% Definiamo una variabile globale affinché possa essere usata da tutte
% le funzioni per il calcolo dei coefficienti e della soluzione
% analitica.
global beta
beta = 10;

% numero di elementi
N = 15;
% grado degli elementi
r = 3;
% posizione dei nodi della griglia: equispaziati
xnod = [0 : 10/(r*N) : 10];

% costruzione della griglia di calcolo
griglia = struttura_griglia(xnod,r);

% definizione della base locale
base = struttura_base(r);

% costruzione del sistema
K = termine_diffusione(griglia,base,'coeff_mu');
C = termine_trasporto(griglia,base,'coeff_a');
% oppure, volendo usare una stabilizzazione:
%C = termine_trasporto_upwind(griglia,base,'coeff_a','coeff_mu');
%C = termine_trasporto_SG(griglia,base,'coeff_a','coeff_mu');
bv = termine_noto(griglia,base,'coeff_f');

% Condizioni al contorno
dati_bordo = struct( ...
  'bc'    , [1 3] , ... % tipo delle condizioni
  'gamma' , [0 3] , ... % parametro gamma
  'r'     , [0 15] ...  % dato
);
[R r] = termine_bordo(griglia,base,dati_bordo);

A = K+C+R;
b = bv+r;

% Condizioni al contorno di prima specie: partizionamento del
% sistemali lineare
u_incognite = [2 : 1 : griglia.dim];
u_note      = [1];

A11 = A(u_incognite,u_incognite);
A12 = A(u_incognite,u_note     );
A21 = A(u_note     ,u_incognite);
A22 = A(u_note     ,u_note     );
b1 = b(u_incognite);
x2 = [0];  % valori noti di u al bordo

% risoluzione
x1 = A11\(b1-A12*x2);

% Soluzione numerica completa: valori calcolati e valori prescritti
x(u_incognite,1) = x1;
x(u_note     ,1) = x2;

% rielaborazione dei risultati e grafica
[x_plot y_plot] = visualizza_soluzione(griglia,base,x,20);
plot(x_plot,sol_esatta(x_plot),'m')

% valutazione dell'errore
[eL2,eH1] = norme_errore(griglia,base,x,'sol_esatta','grad_sol_esatta');
disp(['Errore L2: ',num2str(eL2)])
disp(['Errore H1: ',num2str(eH1)])

% Infine, è possibile calcolare la reazione vincolare in x=0
b2 = A21*x1+A22*x2;
reazione_vincolare = bv(u_note)-b2
