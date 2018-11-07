function [A,b,mod_basis,a_ril,b_ril,yq,wyq]=build_system(size_mb,omega_xl,omega_xr,hx,bc_up,bc_down,dato_uploc,dato_downloc,posizione_dominio,Coeff_forma,Dati,Dati_geometrici,gamma,nd,coupling)

% Numero nodi della formula di quadratura %attenzione, salendo con il numero di modi ?? necessario raffinare la griglia di quadratura
nqnx = 4;
nqny = 64;

% ne numero intervalli e nx numero nodi in direzione X
ne = round((omega_xr-omega_xl)/hx);
nx = ne+1;

% Nodi di intergrazione di Gauss-Legendre
[~, xq, wxq] = gausslegendre( nqnx ); % nqnx non viene cambiato
[~, yq, wyq] = gausslegendre( nqny ); % nqny non viene cambiato

% Mesh FEM in x: nodi equispaziati
meshx = zeros(nx,1);
meshx(1) = omega_xl;
for i=2:nx
    meshx(i) = meshx(i-1)+hx;
end

% Nodi di quadratura lungo tutta la mesh
% Ogni intervallo ha i suoi nodi e i suoi pesi
mesh_xx = zeros( ne*nqnx, 1); % Nodi
mesh_wx = zeros( ne*nqnx, 1); % Pesi

for i=1:ne
    [mesh_xx( (i-1)*nqnx+1 : i*nqnx ), mesh_wx( (i-1)*nqnx+1 : i*nqnx ) ] = quadrature_rule( meshx(i), meshx(i+1), xq, wxq);
end

% Gradi di libert??
% OSS se degree_x=1 size_fb=nx
% altrimenti
% size_fb = (degree_x+1)*ne - (nx-2);
size_fb=nx;

% Base modale per la direzione Y
[mod_basis, mod_basis_y] = new_modal_basis( size_mb, yq, bc_up,bc_down,Coeff_forma,Dati_geometrici.L(0));

% Base degli elementi finiti per la direzione X
[fem_basis, fem_basis_x] = new_fem_basis( 1,meshx, mesh_xx );


% Il numero di elementi non zero della matrice ?? (3*nx-2)*(size_mb+1)^2 (caso P1)

% Matrice di sistema (struttura a blocchi) e termine noto
A = sparse( size_fb*(size_mb), size_fb*(size_mb) );
b = zeros ( size_fb*(size_mb), 1);

[x,YQ] = meshgrid( mesh_xx, yq );
L_computed     = Dati_geometrici.L(x)'; % spessore
a_computed     = Dati_geometrici.a(x)'; % cordinate y inferiore del canale

y = a_computed' + L_computed'.*YQ; % inserisco la dimensionalit??

mu_c = Coeff_forma.mu(x,y)';
if (isa(Coeff_forma.beta1,'function_handle'))
    beta1_c = Coeff_forma.beta1(x,y)';
    beta2_c = Coeff_forma.beta2(x,y)';
else
    beta1_c = Coeff_forma.beta1;
    beta2_c = Coeff_forma.beta2;
end
sigma_c = Coeff_forma.sigma(x,y)';
% 
% nqny=64;
% [~, adhocnodes, adhocw ] = gausslegendre( nqny ); % nqny non viene cambiato
% adhocnodes=adhocnodes/5+0.4;
% adhocw=adhocw/5;
% [xadh,yadh]=meshgrid(mesh_xx,adhocnodes);
% 
%  force_c = Dati.force(xadh,yadh)';
force_c = Dati.force(x,y)';

Computed = struct('mu_c',mu_c,'beta1_c',beta1_c,'beta2_c',beta2_c,'sigma_c',sigma_c,'force_c',force_c,'y',y,'yq',YQ);

% Ciclo di assemblaggio
% Viene considerata ogni coppia di frequenze
% e calcolato il relativo blocco di matrice e di termine noto.
for imb = 1:size_mb
    for kmb = 1:size_mb
        [Amb,bmb,a_ril,b_ril] = assembla_x( imb, kmb, size_fb, 1, ne, nqnx, mesh_wx, wyq, ...
            mod_basis(:,imb), mod_basis_y(:,imb), ...
            mod_basis(:,kmb), mod_basis_y(:,kmb), ...
            fem_basis, fem_basis_x, ...
            L_computed, bc_up, bc_down, ...
            dato_uploc,dato_downloc,posizione_dominio,Coeff_forma,gamma,Computed,nd,coupling);
        
        % Assegnamento del blocco appena assemblato.
        A(1+(imb-1)*size_fb : imb*size_fb , 1+(kmb-1)*size_fb : kmb*size_fb) = Amb;
    end
    b( 1+(imb-1)*size_fb : imb*size_fb ) = bmb;
end