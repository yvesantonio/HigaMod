    function [u]=solver_leg(mesh_x, Coeff_forma, m, force,Map,bc_up,bc_down)
%
% Solver in (0,L)x(-1,1) with legendre polinomial
% Inflow/outflow, homogeneous neumann conditions.
% lateral boundary robin
%

%Setting up a quadrature rule for the y direction.
numberofnodes_y=64;
[trash, node_y, weights_y] = gausslegendre( numberofnodes_y );
% I adjust the nodes and weights to work on (-1,1)
node_y=node_y*2-1;
weights_y=weights_y*2;

%number of nodes
nx=length(mesh_x);
%number of elements
ne=nx-1;

%Setting up a quadrature rule for the x direction.
numberofnodes_x=4;
[trash, node_x, weights_x] = gausslegendre( numberofnodes_x );

quad_nodes_mesh_x = zeros( ne*numberofnodes_x, 1); % Nodi
quad_weights_mesh_x = zeros( ne*numberofnodes_x, 1); % Pesi

for i=1:ne
    [quad_nodes_mesh_x( (i-1)*numberofnodes_x+1 : i*numberofnodes_x ), quad_weights_mesh_x( (i-1)*numberofnodes_x+1 : i*numberofnodes_x ) ] = ...
			 quadrature_rule( mesh_x(i), mesh_x(i+1), node_x, weights_x);
end

[mod_basis, mod_basis_y] = new_modal_basis_legendre( m, node_y,bc_up,bc_down);
[fem_basis, fem_basis_x] = new_fem_basis(1,mesh_x, quad_nodes_mesh_x);

% Il numero di elementi non zero della matrice è (3*nx-2)*(size_mb+1)^2 (caso P1)

% Matrice di sistema (struttura a blocchi) e termine noto
A = sparse( nx*m, nx*m );
b = zeros ( nx*m, 1);

[xhat,yhat] = meshgrid( quad_nodes_mesh_x, node_y);

x = xhat;
y = Map.invhaty(xhat,yhat);

mu_c     = Coeff_forma.mu(x,y)';
beta1_c  = Coeff_forma.beta1(x,y)';
beta2_c  = Coeff_forma.beta2(x,y)';
sigma_c  = Coeff_forma.sigma(x,y)';
force_c  = force(x,y)';
Jac_c    = Map.jac(xhat,yhat)';
D_c      = Map.D(xhat,yhat)';
J_c      = Map.J(xhat,yhat)';
a_down_c = Map.a_down(xhat(1,:))';
a_up_c   = Map.a_up(xhat(1,:))';

jacmu_c=Jac_c.*mu_c;
jacmuD_c=jacmu_c.*D_c;
jacbeta1_c=Jac_c.*beta1_c;
jacmud2j2_c=jacmu_c.*(D_c.^2+J_c.^2);
jacsigma_c=Jac_c.*sigma_c;
jacbeta1Dejacbeta2j_c=jacbeta1_c.*D_c+Jac_c.*beta2_c.*J_c;
jacforce_c=Jac_c.*force_c;

Computed = struct('a_down',a_down_c,'a_up',a_up_c , 'jacforce',jacforce_c ...
    ,'jacmu',jacmu_c,'jacmuD',jacmuD_c,'jacbeta1',jacbeta1_c,'jacmud2j2',jacmud2j2_c,'jacsigma',jacsigma_c,'jacbeta1Dejacbeta2j',jacbeta1Dejacbeta2j_c);

% Ciclo di assemblaggio
% Viene considerata ogni coppia di frequenze
% e calcolato il relativo blocco di matrice e di termine noto.
for imb = 1:m
    for kmb = 1:m
       [Amb,bmb] = assembla_legendre( imb, kmb, nx, quad_weights_mesh_x, weights_y, ...
            mod_basis(:,imb), mod_basis_y(:,imb), ...
            mod_basis(:,kmb), mod_basis_y(:,kmb), ...
            fem_basis, fem_basis_x, ...
            bc_up,bc_down, ...
            Coeff_forma,Computed);
        % Assegnamento del blocco appena assemblato.
        A(1+(imb-1)*nx : imb*nx , 1+(kmb-1)*nx : kmb*nx) = Amb;
    end
    b( 1+(imb-1)*nx : imb*nx ) = bmb;
end
u=A\b;
plot_solution_legendre(m,mesh_x, u,Map,bc_up,bc_down);
export_py_legendre(m,[mesh_x(1),mesh_x(end)],mesh_x(2)-mesh_x(1),u,['m=',num2str(m)],bc_up,bc_down,Map);
axis equal
end
