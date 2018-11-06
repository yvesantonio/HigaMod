%-+--------------------------------------------------------------------
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
function export_py(size_mb,a_ril,b_ril,cutx,hx, u,bc_up,bc_down,Coeff_forma,nome_file,L,a)

%% IMPORT CLASSES
            
            import Core.AssemblerADRHandler
            import Core.BoundaryConditionHandler
            import Core.IntegrateHandler
            import Core.EvaluationHandler
            import Core.BasisHandler
            import Core.SolverHandler
            
%% Possible input problems

numbVerNodes = 100;

if nargin < 11
    L=1;
    a=0;
    M = numbVerNodes;
    warning('a, L and M are not specified, are you sure?')
end

if nargin < 13
   M = numbVerNodes;
    warning('a, L and M are not specified, are you sure?')
end

%% Opens the folder of matlab results

%% Create the Freefem++ simulation folder
    
    for ii = 1:1000
        
        checkFolder = exist(['MatlabOutput',num2str(ii)]);
        
        if (checkFolder ~= 7)
            fileName = ['MatlabOutput',num2str(ii)];
            mkdir(fileName)
            break
        end
    end
    
%% Open the new folder created
    
    cd(fileName);

%% Definition of the number of domains

nd = size(size_mb,2);

%% Definition of the temporary name of the file

temp = 'matlab_result';

%% Computes the number of elements in the X direction

ne = 0;
for i = 1:nd
    ne = ne + round((cutx(i+1)-cutx(i))/hx(i));
end

%% Computes the number of elements in the Y direction

M = numbVerNodes;

%% Adds the file specification to the name

nome_file=strcat(temp,'.out');

%% Creates the file of the Matlab solution

fid =fopen(nome_file,'w');

%% Specifies the number of elements in each direction

% Note: The X direction contains three times the numbers of elements of the
% actual solution because the Freefemm++ simulation is performed using P2
% elements, which implies the use of 2 extra coordinates for each element.

fprintf(fid,'            %u\n',3*ne);
fprintf(fid,'             %u\n',numbVerNodes);

for i=1:nd
    
    %% Number of elements and number of nodes
    
    ne = round((cutx(i+1)-cutx(i))/hx(i));
    nx = ne+1;

    %% Creates the augmented mesh on the X direction to include the nodes
    % used in the P2 Freefem solution
    
    meshx{i} = zeros(3*nx,1);
    meshx{i}(1) = cutx(i);
    meshx{i}(2) = cutx(i)+hx(i)/2;
    meshx{i}(3) = cutx(i)+hx(i);
    
    %% Total number of elements in both directions
    
    numbHorNodes = 3*ne;
    nqny = numbVerNodes;
    
    %% Extremes of the domain along the X direction
    
    omega_xl = cutx(i);
    omega_xr = cutx(i+1);
    
    %% Computes the quadrature nodes to be inserted in the elements along
    % the X direction
  
    obj_gaussLegendre_1 = IntegrateHandler();
    obj_gaussLegendre_1.numbQuadNodes = numbHorNodes;
    [numbHorNodes,xq,wxq] = gaussLegendre(obj_gaussLegendre_1);
    
    %% Computes the quadrature nodes to be inserted in the elements along
    % the Y direction

    obj_gaussLegendre_2 = IntegrateHandler();
    obj_gaussLegendre_2.numbQuadNodes = nqny;
    [nqny,mesh_y,mesh_wy] = gaussLegendre(obj_gaussLegendre_2);
    
    %% Adjust the quadrature nodes along the X direction to consider the
    % actual length of the domain
  
    obj_quadratureRule = IntegrateHandler();

    obj_quadratureRule.leftBoundInterval = omega_xl;
    obj_quadratureRule.rightBoundInterval = omega_xr;
    obj_quadratureRule.inputNodes = xq;
    obj_quadratureRule.inputWeights = wxq;
                
    [mesh_xx, mesh_wx] = quadratureRule(obj_quadratureRule);
    
    %% Add the additional coordinates to the elements of the horizontal mesh
 
    for h=2:ne
        meshx{i}( (h-1)*3 + 1) = meshx{i}( (h-1)*3 );
        meshx{i}( (h-1)*3 + 2) = meshx{i}( (h-1)*3 ) + hx(i)/2;
        meshx{i}( (h-1)*3 + 3) = meshx{i}( (h-1)*3 ) + hx(i);
    end

    meshyref = linspace(0,1,numbVerNodes)';
    meshy = linspace(0,1,numbVerNodes)';
    yvis = length(meshy);
    
    %% Creates the data structure to the solution data to be evaluated on
    % the augmented horizontal mesh

    sol{i}  = zeros(3*nx,yvis);
    solx{i} = zeros(3*nx,yvis);
    soly{i} = zeros(3*nx,yvis);
    
    %% Evaluate the modal basis along the vertical coordinates

    obj_newModalBasis = BasisHandler();
            
    obj_newModalBasis.dimModalBasis = size_mb(i);
    obj_newModalBasis.evalNodesY = meshyref;
    obj_newModalBasis.labelUpBoundCond = bc_up{i};
    obj_newModalBasis.labelDownBoundCond = bc_down{i};
    obj_newModalBasis.coeffForm = Coeff_forma;

    [mb,mby] = newModalBasis(obj_newModalBasis);

    % sol(i,j)=sol(i,j)+(  u(index(i,1)+(imb-1)*nx)*(meshx(index(i,2))-mesh_xx(i))/hx ...
    % index(i,1)==h
    % index(i,2)==k

    %% dx for the first interval, use u(-1/2)=0
    
    h=1;
    for k=1:M
        for imb=1:size_mb(i)
            solx{i}(3*(h-1)+1,k)=solx{i}(3*(h-1) + 1, k )+1/2*(u(h+(imb-1)*nx) - u(h+1+(imb-1)*nx) )/hx*mb(k,imb);
            solx{i}(3*(h-1)+2,k)=solx{i}(3*(h-1) + 2, k )+(u(h+1+(imb-1)*nx)-u(h+(imb-1)*nx))/hx*mb(k,imb);
            solx{i}(3*(h-1)+3,k)=solx{i}(3*(h-1) + 3, k )+(u(h+2+(imb-1)*nx)-u(h+(imb-1)*nx))/(2*hx)*mb(k,imb);
        end
    end



    for h=1:ne
        for k=1:numbVerNodes
            for imb=1:size_mb(i)
                sol{i}(3*(h-1) + 1 , k)=sol{i}(3*(h-1) + 1, k ) + u(h+(imb-1)*nx)*mb(k,imb);
                sol{i}(3*(h-1) + 2 , k)=sol{i}(3*(h-1) + 2, k ) + 1/2*(u(h+(imb-1)*nx) + u(h+1+(imb-1)*nx) ) *mb(k,imb);
                sol{i}(3*(h-1) + 3 , k)=sol{i}(3*(h-1) + 3, k ) + u(h+1+(imb-1)*nx)*mb(k,imb);
                
                if (h>=2 && h<=(ne-1))
                    solx{i}(3*(h-1)+1,k)=solx{i}(3*(h-1) + 1, k ) + 1/2*(u(h+1+(imb-1)*nx) - u(h-1+(imb-1)*nx) )/hx*mb(k,imb);
                    solx{i}(3*(h-1)+2,k)=solx{i}(3*(h-1) + 2, k ) +(u(h+1+(imb-1)*nx)-u(h+(imb-1)*nx))/hx*mb(k,imb);
                    solx{i}(3*(h-1)+3,k)=solx{i}(3*(h-1) + 3, k ) +(u(h+2+(imb-1)*nx)-u(h+(imb-1)*nx))/(2*hx)*mb(k,imb);
                end 
            
                soly{i}(3*(h-1) + 1 , k)=soly{i}(3*(h-1) + 1, k ) + u(h+(imb-1)*nx)*mby(k,imb);
                soly{i}(3*(h-1) + 2 , k)=soly{i}(3*(h-1) + 2, k ) + 1/2*(u(h+(imb-1)*nx) + u(h+1+(imb-1)*nx) ) *mby(k,imb);
                soly{i}(3*(h-1) + 3 , k)=soly{i}(3*(h-1) + 3, k ) + u(h+1+(imb-1)*nx)*mby(k,imb);               
            end
            
            sol{i}(3*(h-1) + 1,k)=sol{i}(3*(h-1)+1,k)+ a_ril(i)*meshy(k)+b_ril(i);
            sol{i}(3*(h-1) + 2,k)=sol{i}(3*(h-1)+2,k)+ a_ril(i)*meshy(k)+b_ril(i);
            sol{i}(3*(h-1) + 3,k)=sol{i}(3*(h-1)+3,k)+ a_ril(i)*meshy(k)+b_ril(i);
           
            soly{i}(3*(h-1) + 1 , k)=soly{i}(3*(h-1) + 1, k ) + a_ril(i)*meshy(k)+b_ril(i);
            soly{i}(3*(h-1) + 2 , k)=soly{i}(3*(h-1) + 2, k ) + a_ril(i)*meshy(k)+b_ril(i);
            soly{i}(3*(h-1) + 3 , k)=soly{i}(3*(h-1) + 3, k ) + a_ril(i)*meshy(k)+b_ril(i);
        end
    end

    %% dx for the last interval, use u(6+0.5)=0
    
    for k=1:numbVerNodes
        for imb=1:size_mb(i)
            solx{i}(3*(h-1)+1,k)=solx{i}(3*(h-1) + 1, k ) + 1/2*(u(h+1+(imb-1)*nx) - u(h-1+(imb-1)*nx) )/hx*mb(k,imb);
            solx{i}(3*(h-1)+2,k)=solx{i}(3*(h-1) + 2, k ) +(u(h+1+(imb-1)*nx)-u(h+(imb-1)*nx))/hx*mb(k,imb);
            solx{i}(3*(h-1)+3,k)=solx{i}(3*(h-1) + 3, k ) +(-u(h+(imb-1)*nx))/(2*hx)*mb(k,imb);
        end
    end
    
    %% Writes the computed solution in the simulation file

    for k=1:ne
        for j=1:M
            fprintf(fid,'       %1.7E       %1.7E    %1.7E   %1.7E   %1.7E\n',meshx{i}( (k-1)*3 + 1),meshy(j),sol{i}((k-1)*3 + 1,j),solx{i}((k-1)*3 + 1,j),soly{i}((k-1)*3 + 1,j)   );
        end
        for j=1:M
            fprintf(fid,'       %1.7E       %1.7E    %1.7E  %1.7E  %1.7E\n',meshx{i}( (k-1)*3 + 2),meshy(j),sol{i}((k-1)*3 + 2,j),solx{i}((k-1)*3 + 2,j),soly{i}((k-1)*3 + 2,j) );
        end
        for j=1:M
            fprintf(fid,'       %1.7E       %1.7E     %1.7E  %1.7E  %1.7E\n',meshx{i}( (k-1)*3 + 3),meshy(j),sol{i}((k-1)*3 + 3,j),solx{i}((k-1)*3 + 3,j),soly{i}((k-1)*3 + 3,j) );
        end
    end
    
    %% Closes the simulation files
    
    fclose(fid);
    
    %% Creates the name of the file to collect the quadrature nodes
    
    temp=['quadrature_node'];
    nome_file=strcat(temp,'.out');
    
    %% Creates the file that containing the quadrature nodes
    
    fid =fopen(nome_file,'w');
    fprintf(fid,'            %u\n',numbHorNodes);
    fprintf(fid,'            %u\n',nqny);
    num2str(size_mb)
    
    %% Writes the quadrature nodes and weights for the original mesh and the
    % augmented mesh in the X direction
    
    for k=1:numbHorNodes
       fprintf(fid,'   %1.7E  %1.7E   %1.7E    %1.7E  \n', xq(k), wxq(k), mesh_xx(k) , mesh_wx(k));
    end

    %% Writes the quadrature nodes and weights in the Y direction
    
    for k=1:numbVerNodes
       fprintf(fid,'   %1.7E  %1.7E       \n', mesh_y(k),mesh_wy(k));
    end
    
    %% Writes the number of nodes and used in each direction
    
    fprintf(fid,'   %1.7E  %1.7E     \n',numbHorNodes,numbVerNodes);
    
end

%% Closes the file containing the quadratures nodes

fclose(fid);

%% Exit the Matlab results folder

cd ..

return
