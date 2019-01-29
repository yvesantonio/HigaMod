%----------------------------------------------------------------------
% HiModlab is a general purpose Hierarchical model reduction library.
% Copyright (C) 2006-2017  by the HiMod authors (see authors.txt).
%

% This file is part of the HiModlab library.
%
% The HiModlab library is free software: you can use it, redistribute
% it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either
% version 2 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Acknowledgement and Disclaimer
%
% The development of HiModLab was supported by the US National Science
% Foundation (NSF Project, DMS 1419060 ). Any opinions, findings and conclusions
% or recomendations expressed in the source code and documentation are
% those of the authors and do not necessarily reflect the views of the
% National Science Foundation (NSF).
%
%
%
%-+--------------------------------------------------------------------
%% HI-MOD MONODOMAIN
% The following script allows the solution of an Advection - Diffusion -
% Reaction differential problem in 2D using the HIGAMod solution.

    clear all
    close all

    disp('******************************************')
    disp('*           HIGAMod Simulation           *');
    disp('******************************************');

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

    %% Simulation case   
    caso  = 2;      % Analysed Case
    
    %-------------------------------------------------------------------------%
    % Note;
    % The exciting force acting uppon the system is different depending on the
    % case we are analysing. In the current version of the code we are
    % considering the following exciting forces:
    %
    %  (1) :  Straight centerline;
    %  (2) :  Parabolic centerline;
    %  (3) :  Cubic centerline;
    %-------------------------------------------------------------------------%

    switch caso
        
    case {1,2,3,4,5,6,7,8,9,10}
        minHor     = 0;
        maxHor     = 1;
        minVer     = 0;
        maxVer     = 1;
    end

    %-------------------------------------------------------------------------%
    % NOTE:
    % If one wants to change the values of min_x and max_x, the procedure can
    % be performed inside the class 'AssemblerIGA'.
    %-------------------------------------------------------------------------%

    %% Discretization parameters
    
    discStruct = [];
    
    modesVelc = 10;
    modesPres = 5;
    
    % Step used to generate the knot vector
    
    discStruct.stepHorMesh     = 0.1;
    
    % Number of knots/elements in the isogeometric space
    
    discStruct.numbElements    = 1/discStruct.stepHorMesh;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    discStruct.numbModesP       = modesPres;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE X COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of transverse modes
    
    discStruct.numbModesUx       = modesVelc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE Y COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of transverse modes
    
    discStruct.numbModesUy       = modesVelc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISCRETIZATION PARAMETERS FOR THE Z COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of transverse modes
    
    discStruct.numbModesUz       = modesVelc;
    
    %% Isogeometric basis properties
    
    igaBasisStruct = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Polynomial Degree of the B-Spline Base
    
    igaBasisStruct.degreeSplineBasisP    = 1;
    
    % Continuity of the Base 'C^(p-k)'
    
    igaBasisStruct.continuityParameterP  = 0;
    
    % Number of control points
    
    igaBasisStruct.numbControlPtsP = discStruct.numbElements * ...
                                     (igaBasisStruct.degreeSplineBasisP - igaBasisStruct.continuityParameterP) + ...
                                     1 + igaBasisStruct.continuityParameterP;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE X COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Polynomial Degree of the B-Spline Base
    
    igaBasisStruct.degreeSplineBasisUx    = 2;
    
    % Continuity of the Base 'C^(p-k)'
    
    igaBasisStruct.continuityParameterUx  = 1;
    
    % Number of control points
    
    igaBasisStruct.numbControlPtsUx = discStruct.numbElements * ...
                                      (igaBasisStruct.degreeSplineBasisUx - igaBasisStruct.continuityParameterUx) + ...
                                      1 + igaBasisStruct.continuityParameterUx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE Y COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Polynomial Degree of the B-Spline Base
    
    igaBasisStruct.degreeSplineBasisUy    = 2;
    
    % Continuity of the Base 'C^(p-k)'
    
    igaBasisStruct.continuityParameterUy  = 1;
    
    % Number of control points
    
    igaBasisStruct.numbControlPtsUy = discStruct.numbElements * ...
                                      (igaBasisStruct.degreeSplineBasisUy - igaBasisStruct.continuityParameterUy) + ...
                                      1 + igaBasisStruct.continuityParameterUy;
                                  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IGA BASIS FOR THE Z COMPONENT OF THE VELOCITY %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Polynomial Degree of the B-Spline Base
    
    igaBasisStruct.degreeSplineBasisUz    = 2;
    
    % Continuity of the Base 'C^(p-k)'
    
    igaBasisStruct.continuityParameterUz  = 1;
    
    % Number of control points
    
    igaBasisStruct.numbControlPtsUz = discStruct.numbElements * ...
                                      (igaBasisStruct.degreeSplineBasisUz - igaBasisStruct.continuityParameterUz) + ...
                                      1 + igaBasisStruct.continuityParameterUz;
    
    %% Boundary conditions
    
    boundCondStruct = [];
    
    %%%%%%%%%%%%%
    % PARAMTERS %
    %%%%%%%%%%%%%
    
    % INFLOW PRESSURE
    
    Pinf = 10;
    
    % OUTFLOW PRESSURE
    
    Pout = 0;
    
    % PRESSURE DROP
    
    dP = Pinf - Pout;
    
    % DYNAMIC VISCOSITY
    
    nu = 0.1;
        
    % LENGTH OF THE CYLINDER
    
    L = 5;
    
    % RADIUS OF THE CYLINDER
    
    R = 0.5;
    
    % PARABOLIC PROFILE PARAMETERS
    
    a0 = 0;
    a1 = (dP * R) / (2 * nu * L);
    a2 = (-1 * dP) / (4 * nu * L);
    
    Umax = a1^2 / 4 * a2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR THE PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    boundCondStruct.bc_up_tag_P    = 'neu';
    boundCondStruct.bc_down_tag_P  = 'neu';
    boundCondStruct.bc_inf_tag_P   = 'neu';
    boundCondStruct.bc_out_tag_P   = 'neu';
    boundCondStruct.bc_up_data_P   = 0;
    boundCondStruct.bc_down_data_P = 0;
    boundCondStruct.bc_inf_data_P  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;
    boundCondStruct.bc_out_data_P  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR Ux %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    boundCondStruct.bc_up_tag_Ux    = 'dir';
    boundCondStruct.bc_down_tag_Ux  = 'dir';
    boundCondStruct.bc_inf_tag_Ux   = 'dir';
    boundCondStruct.bc_out_tag_Ux   = 'neu';
    boundCondStruct.bc_up_data_Ux   = 0;
    boundCondStruct.bc_down_data_Ux = 0;
    boundCondStruct.bc_inf_data_Ux  = @(rho,eta) a0 + a1 * rho + a2 * rho.^2;
    boundCondStruct.bc_out_data_Ux  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR Uy %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    boundCondStruct.bc_up_tag_Uy    = 'dir';
    boundCondStruct.bc_down_tag_Uy  = 'dir';
    boundCondStruct.bc_inf_tag_Uy   = 'dir';
    boundCondStruct.bc_out_tag_Uy   = 'dir';
    boundCondStruct.bc_up_data_Uy   = 0;
    boundCondStruct.bc_down_data_Uy = 0;
    boundCondStruct.bc_inf_data_Uy  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;
    boundCondStruct.bc_out_data_Uy  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INFLOW, OUTFLOW AND LATERAL BOUNDARY CONDITIONS FOR Uz %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    boundCondStruct.bc_up_tag_Uz    = 'dir';
    boundCondStruct.bc_down_tag_Uz  = 'dir';
    boundCondStruct.bc_inf_tag_Uz   = 'dir';
    boundCondStruct.bc_out_tag_Uz   = 'dir';
    boundCondStruct.bc_up_data_Uz   = 0;
    boundCondStruct.bc_down_data_Uz = 0;
    boundCondStruct.bc_inf_data_Uz  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;
    boundCondStruct.bc_out_data_Uz  = @(rho,eta) 0 + 0 * rho + 0 * rho.^2;
    
    %% Physical domain
    %---------------------------------------------------------------------%
    % Note: Complete domain is defined using the nurbs functions, not only
    % the centreline. This way, we can automatically compute the map, its
    % first and second derivatives and the jacobian of the transformation
    % from the pysical domain to the reference domain, where the reduction
    % procedure is defined.
    %---------------------------------------------------------------------%

    switch caso
        
    case {1,4,5,6,7,8,9,10}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Artery from PATIENT 1 %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        currentFolder = pwd;
        
        folder = [pwd,'/Geometry/3D_Patient_1/'];
        cd(folder);
        
        filename1 = '1Patient3DGeo.mat';
        filename2 = '1Patient3DMap.mat';
        
        load(filename1);
        load(filename2);
        
        cd ..
        cd ..
        
        Vol = geo;
        
        % figure
        % nrbplot(Vol,[15 15 500]);
        
    case{2}
        
        %%%%%%%%%%%%
        % Cylinder %
        %%%%%%%%%%%%
        
        currentFolder = pwd;
        
        folder = [pwd,'/Geometry/3D_Cylinder/'];
        cd(folder);
        
        filename1 = 'Cylinder3DGeo.mat';
        filename2 = 'Cylinder3DMap.mat';
        
        load(filename1);
        load(filename2);
        
        cd ..
        cd ..
        
        Vol = geo;
        
        % figure
        % nrbplot(Vol,[10 10 10]);
        
    case {3}
        
        %%%%%%%%
        % Slab %
        %%%%%%%%
        
        currentFolder = pwd;
        
        folder = [pwd,'/Geometry/3D_Slab/'];
        cd(folder);
        
        filename1 = 'Slab3DGeo.mat';
        filename2 = 'Slab3DMap.mat';
        
        load(filename1);
        load(filename2);
        
        cd ..
        cd ..
        
        Vol = geo;

    end
    
    %% Quadrature propeties
    %---------------------------------------------------------------------%
    % Specifies the number of quadrature nodes to be used on the horizontal
    % and vertical direction to computes the integrals inside the Build
    % class.
    %---------------------------------------------------------------------%
    
    quadProperties = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pressure quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbHorNodesP = 16;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pressure quadrature nodes along Y direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbVerNodesP = 64;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected X quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbHorNodesUx = 16;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected X quadrature nodes along Y direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbVerNodesUx = 64;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected Y quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbHorNodesUy = 16;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected Y quadrature nodes along Y direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbVerNodesUy = 64;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected Z quadrature nodes along X direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbHorNodesUz = 16;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity field projected Z quadrature nodes along Z direction %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    quadProperties.numbVerNodesUz = 64;
    
    %% Problem parameters

    probParameters = [];
    
    switch caso
    case {1,2,3,4,5,6,7,8,9,10}
        
        % Pressure compensation
        
        probParameters.delta    = @(x,y,z) ( 0.0 + 0*x + 0*y + 0*z );
        
        % Artificial reaction
        
        probParameters.sigma    = @(x,y,z) ( 0.0 + 0*x + 0*y + 0*z );
        
        % Fluid kinetic viscosity
        
        probParameters.nu    = @(x,y,z) (  nu + 0*x + 0*y + 0*z );
        
        % Forcing term acting on the fluid

        probParameters.force_x = @(x,y,z) (  0.00 + 0*x + 0*y + 0*z );
        probParameters.force_y = @(x,y,z) (  0.00 + 0*x + 0*y + 0*z );
        probParameters.force_z = @(x,y,z) (  0.00 + 0*x + 0*y + 0*z );
        
        % Robin coefficient for the modal basis
        
        probParameters.coeffrobin = 1;
        
    end
    
    %% Solver
    
    % Definition of the Object of the EvaluationHandler Class

    import Core.SolverHandler

    obj_solverIGA = SolverHandler();

    % Properties Assignment

    obj_solverIGA.discStruct = discStruct;
    obj_solverIGA.boundCondStruct = boundCondStruct;
    obj_solverIGA.igaBasisStruct = igaBasisStruct;
    obj_solverIGA.probParameters = probParameters;
    obj_solverIGA.geometricInfo = geometricInfo;
    obj_solverIGA.quadProperties = quadProperties;

    tic;
    [errStruct] = solverIGAScatterStokes3D(obj_solverIGA);
    toc;