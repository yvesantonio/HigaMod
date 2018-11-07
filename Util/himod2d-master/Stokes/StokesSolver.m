function sol=StokesSolver(mesh,quad,Data,out,MapFun,steady)

% structuring the mesh
ne=mesh.nb_elem;
nx=ne+1;
mesh_x=linspace(mesh.I(1),mesh.I(2),nx);
hx=mesh_x(2)-mesh_x(1);

%Setting up a quadrature rule for the y direction.
[~, node_y, weights_y] = gausslegendre( quad.ny );
% I adjust the nodes and weights to work on (-1,1)
node_yhat=node_y*2-1;
weights_yhat=weights_y*2;

%Setting up a quadrature rule for the x direction.
[~, node_x, weights_x] = gausslegendre( quad.nx );
quad_nodes_mesh_x = zeros( ne*quad.nx, 1);
quad_weights_mesh_x = zeros( ne*quad.nx, 1);
for i=1:ne
    [ ...
        quad_nodes_mesh_x( (i-1)*quad.nx+1 : i*quad.nx ), ...
        quad_weights_mesh_x( (i-1)*quad.nx+1 : i*quad.nx ) ...
        ] = quadrature_rule( mesh_x(i), mesh_x(i+1), node_x, weights_x);
end
[Xhat,Yhat]=meshgrid(quad_nodes_mesh_x,node_yhat);

%some data
fedof=nx+ne;
ndof=Data.m*(fedof);
mp=Data.m-2;
fedofP=nx;
ndofP=nx*mp;

%Calcola le basi FEM
[ FemBasisU , DFemBasisU ]=fembasis(2,mesh_x(1:2),quad_nodes_mesh_x(1:quad.nx));
FemBasisP = fembasis(1,mesh_x(1:2),quad_nodes_mesh_x(1:quad.nx));

%Calcola le basi Modali
[ModalBasisX,DModalBasisX]=modalbasis_legendre(Data.m, node_yhat, 'dir', 'dir');
[ModalBasisY,DModalBasisY]=modalbasis_legendre(Data.m, node_yhat, 'dir', 'dir');
ModalBasisP=modalbasis_legendre(mp, node_yhat, 'rob', 'rob');

t=0;
uold_x=zeros(ndof,1);
uold_y=zeros(ndof,1);
if(steady)
    Data.ntimestep=1;
end
for t_iter=1:Data.ntimestep
    
    t = t + Data.dt;
    display([ 'Starting ',num2str(t_iter),...
        ' time step out of ', num2str(Data.ntimestep),...
        ' t = ', num2str(t)])
    
    %Valuta la mappaf
    Jac=MapFun.jac(t,Xhat,Yhat)';
    D=MapFun.D(t,Xhat,Yhat)';
    J=MapFun.J(t,Xhat,Yhat)';
    [Xinout,Y]=meshgrid(mesh.I,node_yhat);
    Jin_out=MapFun.jac(t,Xinout,Y)';
    Map=struct('Jac',Jac,'D',D,'J',J);
    Y=MapFun.invhaty(t,Xhat,Yhat);
    W=MapFun.w(t,Xhat,Yhat)';
    %display(['max w  =', num2str(max(max(W)))])
    
    %Valuta la velocità
    Uoldx=EvalHimod(Xhat,uold_x,ModalBasisX,FemBasisU,fedof,quad.nx);
    Uoldy=EvalHimod(Xhat,uold_y,ModalBasisY,FemBasisU,fedof,quad.nx);
    %valuta la forzante
    Fx=Data.forcex(Xhat,Y,t)' + Uoldx/Data.dt;
    Fy=Data.forcey(Xhat,Y,t)' + Uoldy/Data.dt;
    
    
    % Assemblaggio%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ax = sparse( ndof, ndof );
    Ay = sparse( ndof, ndof );
    Ayx = sparse( ndof, ndof );
    Axy = sparse( ndof, ndof );
    
    Bx = sparse( ndofP,ndof );
    By = sparse( ndofP,ndof );
    
    bx = zeros( ndof ,1);
    by = zeros( ndof ,1);
    
    for k = 1:Data.m
        
        MbKX=ModalBasisX(:,k);
        DMbKX=DModalBasisX(:,k);
        
        MbKY=ModalBasisY(:,k);
        DMbKY=DModalBasisY(:,k);
        
        for r = 1:Data.m
            
            MbRX=ModalBasisX(:,r);
            DMbRX=DModalBasisX(:,r);
            
            MbRY=ModalBasisY(:,r);
            DMbRY=DModalBasisY(:,r);
            if(~steady)
                [ Ax(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof), ...
                    bx(1+(r-1)*fedof : r*fedof ) ]= ...
                    AssembleVel1DX(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                    quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',1/Data.dt,Data.nu,Fx,Uoldx,Uoldy-W);
            else
                [ Ax(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof), ...
                    bx(1+(r-1)*fedof : r*fedof ) ]= ...
                    AssembleVel1DX(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                    quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',0,Data.nu,Fx,Uoldx,Uoldy-W);
            end
            Axy(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof)=...
                AssembleVelXY(Map,weights_yhat,...
                DMbKX,MbRY,DMbRY,...
                FemBasisU',DFemBasisU',ne,...
                quad_weights_mesh_x,quad.nx,Data.nu);
            
            Ayx(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof)=...
                AssembleVelYX(Map,weights_yhat,...
                MbKY,DMbKY,DMbRX,...
                FemBasisU',DFemBasisU',ne,...
                quad_weights_mesh_x,quad.nx,Data.nu);
            if(~steady)
            [ Ay(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof), ...
                by(1+(r-1)*fedof : r*fedof ) ]= ...
                AssembleVel1DY(Map,weights_yhat, MbKY, MbRY,DMbKY,DMbRY,ne,...
                quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',1/Data.dt,Data.nu,Fy,Uoldx,Uoldy-W);
            else
                [ Ay(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof), ...
                by(1+(r-1)*fedof : r*fedof ) ]= ...
                AssembleVel1DY(Map,weights_yhat, MbKY, MbRY,DMbKY,DMbRY,ne,...
                quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',0,Data.nu,Fy,Uoldx,Uoldy-W);
           
            end
        end
    end
    for k=1:Data.m
        MbKX=ModalBasisX(:,k);
        DMbKX=DModalBasisX(:,k);
        
        DMbKY=DModalBasisY(:,k);
        
        for s=1:mp
            MbSP=ModalBasisP(:,s);
            Bx(1+(s-1)*fedofP:s*fedofP, 1+(k-1)*fedof : k*fedof) = ...
                AssembleBX(Map,weights_yhat,MbKX,MbSP,DMbKX,ne,quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',FemBasisP');
            By(1+(s-1)*fedofP:s*fedofP, 1+(k-1)*fedof : k*fedof) = ...
                AssembleBY(Map,weights_yhat,MbSP,DMbKY,ne,quad_weights_mesh_x,quad.nx,FemBasisU',FemBasisP');
        end
    end
    [Ax,Ay,bx,by]=BChandler(Ax,Ay,bx,by,Data.BC,Data.m,Data.m,fedof,weights_yhat,ModalBasisX,Jin_out,t);
    % % vecchia maniera di farlo.
    %     %BC per u in inflow
    %     if(Data.essential_inflow)
    %         [Ax,bx]=DirichletInflow(Ax,bx,Data.inflow_dirichlet,Data.m,fedof);
    %     else
    %         bx=NeumannInflow(weights_yhat,bx,Data.Pinf(t),ModalBasisX,fedof,Jin_out');
    %     end
    %
    %     % Imposing uy=0 in outflow
    %     [Ay,by]=DirichletInflow(Ay,by,zeros(Data.m,1),Data.m,fedof); %uy=0 all'inflow
    %     [Ay,by]=DirichletOutflow(Ay,by,zeros(Data.m,1),Data.m,fedof); %uy=0 all'outflow
    
    % FINE ASSEMBLAGGIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=[Ax,Ayx,Bx';Axy,Ay,By';Bx,By,zeros(ndofP,ndofP)];
%     spy(A)
    
    % Soluzione sistema lineare
    sol=A\[bx;by;zeros(ndofP,1)];
    
    uold_x=sol(1:ndof);
    uold_y=sol(ndof:2*ndof);
    
    % Post Processing
    if (t_iter==Data.ntimestep)
        out.ended=true;
    end
    FEnodes=mesh_x(1):hx/2:mesh_x(end);
    [Xnodeshat,Ynodeshat]=meshgrid(FEnodes,node_yhat);
    Yout=MapFun.invhaty(t,Xnodeshat,Ynodeshat);
    plot_stokes(sol,mesh_x,Yout,ModalBasisX,ModalBasisY,ModalBasisP,out,t,0*FEnodes,0*FEnodes,Data.m,Data.m,mp);
end
