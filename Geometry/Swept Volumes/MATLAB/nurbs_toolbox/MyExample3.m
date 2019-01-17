%% DEFINITION OF THE PROBLEM DATA %%
% Note: Here we define the geometry to be approximated, the force, boundary
% conditions and lifting operators.
clc
clear all
close all
clear problem_data  

% PHYSICAL DOMAIN
% Note: Definition of the geometry file to be used in the representation of
% the domain

% HELIX
% coefs =[ 6.0  0.0  6.0  1; 
%         -5.5  0.5  5.5  1; 
%         -5.0  1.0 -5.0  1; 
%          4.5  1.5 -4.5  1; 
%          4.0  2.0  4.0  1; 
%         -3.5  2.5  3.5  1; 
%         -3.0  3.0 -3.0  1; 
%          2.5  3.5 -2.5  1; 
%          2.0  4.0  2.0  1; 
%         -1.5  4.5  1.5  1; 
%         -1.0  5.0 -1.0  1; 
%          0.5  5.5 -0.5  1; 
%          0.0  6.0  0.0  1]'; 
% knots = [0 0 0 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1 1 1]; 
%  
% crv = nrbmak(coefs,knots); 
% nrbplot(crv,100); 
% title('3D helical curve.'); 

% CIRCUNFERENCE ARC
r = 1;
problem_data.geo_name = nrbcirc(r,[],deg2rad(0),deg2rad(180));

% STRAIGHT LINE
%problem_data.geo_name = nrbline ([0 0], [1 1]);

% BOUNNDARY CONDITIONS
% Note: Definition of the boundary conditions in each boundary of the
% problem.

problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2];

% PROBLEM PARAMETERS
% Note: Diffusion coefficient, reaction coefficient and advective field
% used in the definition of the ADR problem.

problem_data.c_diff  = @(x,y) 1 + 0 * x + 0 * y;
problem_data.c_reac  = @(x,y) 1 + 0 * x + 0 * y;

% FORCING TERM

mu = @(x,y) problem_data.c_diff(x,y);
sigma = @(x,y) problem_data.c_reac(x,y);

problem_data.f = @(x,y) (mu(x,y) + sigma(x,y)).* y ./ r;

% LIFTING OPERATORS

problem_data.g = @(x,y, ind) 0;
problem_data.h = @(x,y, ind) 0;

% EXACT SOLUTION
% Note: Note required to solve the problem, but allows the direct
% computational of the approximation error

problem_data.uex     = @(x,y) y ./ r;
problem_data.graduex = @(x,y) x ./ r;

%% CHOICE OF THE DISCRETIZATION PARAMETERS %%
% Note: definition of the parameters used in the isogeometric analysis.
% Power of the basis functions, knot vector, continuity properties in each
% patch, etc.

clear method_data

method_data.degree     = 3;     % Degree of the splines
method_data.regularity = 2;     % Regularity of the splines
method_data.nsub       = 10;    % Number of subdivisions
method_data.nquad      = 4;     % Points for the Gaussian quadrature rule

%% STIFFNESS MATRIX AND RHS ASSEMBLER %%

% EXTRACT FIELDS FROM DATA STRUCTURES 

data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% DOMAIN PROPERTIES
% Note: Construct the geometry properties from the geometric file. We use
% the reference geometry and insert the required number of subdivisions in
% each patch to refine the approximation of the solution.

geometry  = geo_load (geo_name);

[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

% CONSTRUCT THE MESH
% Note: Use the domain property to create the 1D mesh and include the
% quadrature nodes necessary to compute the integration of the elements of
% the stiffness matrix.

rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% CONSTRUCT THE FUNCTINAL SPACE
% Note: Construct the functional space with the basis functions to be used
% in the approximation of the differential problem.

space    = sp_bspline (knots, degree, msh);

% ASSEMBLE THE STIFFNESS MATRIX
% Note: Use the precomputed operator funstions to compute the elements of
% the stiffness matrix.

stiff_mat = op_gradu_gradv_tp (space, space, msh, c_diff) + ...
            op_u_v_tp (space, space, msh, c_reac);
rhs       = op_f_v_tp (space, msh, f);

% NEUMANN BOUNDARY CONDITIONS
% Note: Modify the orinial stiffness matrix and right-handside term to
% include the Neumann boundary condtions, if any.

for iside = nmnn_sides
  if (msh.ndim > 1)
% Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
    gside = @(varargin) g(varargin{:},iside);
    dofs = space.boundary(iside).dofs;
    rhs(dofs) = rhs(dofs) + op_f_v_tp (space.boundary(iside), msh.boundary(iside), gside);
  else
    if (iside == 1)
      x = msh.breaks{1}(1);
    else
      x = msh.breaks{1}(end);
    end
    sp_side = space.boundary(iside);
    rhs(sp_side.dofs) = rhs(sp_side.dofs) + g(x,iside);
  end
end

% DIRICHLET BOUNDARY CONDITIONS
% Note: Modify the orinial stiffness matrix and right-handside term to
% include the Dirichlet boundary condtions, if any.

u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% SOLVE THE LINEAR SYSTEM

u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

%% POST PROCESSING OF RESULTS
% Note: Export the solution file in VTK format readable by Paraview and
% compute the error using the exact solution.

% PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

vtk_pts = {linspace(0, 1, 1000)};

[eu, F] = sp_eval (u, space, geometry, 10);
% plot (F, eu, F, problem_data.uex(F))
plot3(F(1,:),F(2,:),eu,'o'); hold on;
plot3(F(1,:),F(2,:),zeros(size(eu))); hold on;
plot3(F(1,:),F(2,:),problem_data.uex(F(1,:),F(2,:)));
legend ('Numerical solution','Domain', 'Exact solution'); axis tight

% EXPORT PARAVIEW FILE

output_file = 'MyExample3.vts';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% ERROR IN L2 AND H1 NORM

% [error_h1, error_l2] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);