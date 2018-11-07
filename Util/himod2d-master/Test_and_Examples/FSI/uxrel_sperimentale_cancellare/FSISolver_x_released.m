function A=FSISolver_x_released(mesh,quad,Data,out)

%structuring the mesh
ne=mesh.nb_elem;
nx=ne+1;
mesh_x=linspace(mesh.I(1),mesh.I(2),nx);
% hx=mesh_x(2)-mesh_x(1);

%Setting up a quadrature rule for the y direction.
[~, node_y, weights_y] = gausslegendre( quad.ny );
% I adjust the nodes and weights to work on (-1,1)
node_yhat=node_y*2-1;
weights_yhat=weights_y*2;

%Setting up a quadrature rule for the x direction.
[~, node_x, weights_x] = gausslegendre( quad.nx);
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
fedofP=nx;
fedof=nx+ne;
mx=Data.m-2;
my=Data.m;
mp=Data.m-2;
ndofy=my*fedof;
ndofx=mx*fedof;
ndofP=mp*fedofP;

%Calcola le basi FEM
[ FemBasisU , DFemBasisU ]=fembasis(2,mesh_x(1:2),quad_nodes_mesh_x(1:quad.nx));
FemBasisP = fembasis(1,mesh_x(1:2),quad_nodes_mesh_x(1:quad.nx));

%Calcola le basi Modali
[ModalBasisX,DModalBasisX]=modalbasis_legendre(mx, node_yhat, 'dir', 'dir');
[ModalBasisY,DModalBasisY]=modalbasis_legendre(my, node_yhat, 'rob', 'rob');
[MBrobY,~]=modalbasis_legendre(my, [1,-1], 'rob', 'rob');
ModalBasisP=modalbasis_legendre(mp, node_yhat, 'rob', 'rob');
[ModalBasisE,DModalBasisE]=modalbasis_legendre(2, node_yhat,'rob','rob');

ModalBasisX=[ModalBasisE,ModalBasisX];
DModalBasisX=[DModalBasisE,DModalBasisX];

%Inizializzazioni
t=0;
PF_max=30;

%Velocity
uold_x=zeros(ndofx+2,1);
uold_y=zeros(ndofy,1);
u_kmeno1_x=uold_x;
u_kmeno1_y=uold_y;

%spostamento
etaup_old=zeros(fedof,1);
etadown_old=zeros(fedof,1);
etaup_k=etaup_old;
etadown_k=etadown_old;
%Mappae e Jacobiani
Jac=ones(ne*quad.nx,quad.ny);
D=zeros(ne*quad.nx,quad.ny);
J=ones(ne*quad.nx,quad.ny);
Jin_out=ones(quad.ny,2);
Aup=ones(ne*quad.nx,1);
Adown=ones(ne*quad.nx,1);
Map=struct('Jac',Jac,'D',D,'J',J,'Aup',Aup,'Adown',Adown);
Y=Yhat;
Yold=Y;
W=zeros(ne*quad.nx,quad.ny);

% Ciclo in tempo
for t_iter=1:Data.ntimestep
    
    test=1;
    PF_iter=0;
    t = t + Data.dt;
    display([ 'Starting ',num2str(t_iter),...
        ' time step out of ', num2str(Data.ntimestep),...
        ' t = ', num2str(t)])
    
    %Valuta la velocità old
    Uoldx=EvalHimodREL(Xhat,uold_x,ModalBasisX,FemBasisU,fedof,quad.nx);
    Uoldy=EvalHimod(Xhat,uold_y,ModalBasisY,FemBasisU,fedof,quad.nx);
    
    % Iterazioni di punto fisso
    while(test>Data.toll && PF_iter<PF_max)
        
        Ax = sparse( ndofx+2, ndofx+2 );
        Ay = sparse( ndofy, ndofy );
        
        Axy = sparse( ndofy, ndofx +2);
        Ayx = sparse( ndofx+2, ndofy );
        
        Bx = sparse( ndofP,ndofx +2);
        By = sparse( ndofP,ndofy );
        
        bx = zeros( ndofx+2 ,1);
        by = zeros( ndofy ,1);
        
        PF_iter=PF_iter+1;
        %valuta la forzante
        Fx=Data.forcex(Xhat,Y,t)' + Uoldx/Data.dt;
        Fy=Data.forcey(Xhat,Y,t)' + Uoldy/Data.dt;
        
        %valuta il termine non lineare per il punto fisso
        U_kmeno1_x=EvalHimodREL(Xhat,u_kmeno1_x,ModalBasisX,FemBasisU,fedof,quad.nx);
        U_kmeno1_y=EvalHimod(Xhat,u_kmeno1_y,ModalBasisY,FemBasisU,fedof,quad.nx);
        
        %Valuta lo spostamento della struttura
        [Etaup_old,Etadown_old]=EvalDisplacement(etaup_old,etadown_old,FemBasisU,Xhat(1,:),quad.nx);
        %Assemblaggio
        for k = 1:mx+2
            
            MbKX=ModalBasisX(:,k);
            DMbKX=DModalBasisX(:,k);
            
            for r = 1:mx+2
                MbRX=ModalBasisX(:,r);
                DMbRX=DModalBasisX(:,r);
                if(k>2&&r>2)
                    [ Ax(1+(r-1-2)*fedof +2: (r-2)*fedof+2 , 1+(k-1-2)*fedof +2: (k-2)*fedof+2), ...
                        bx(1+(r-1-2)*fedof+2 : (r-2)*fedof +2) ]= ...
                        AssembleVel1DXREL(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                        quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',Data.dt,Data.nu,Fx,U_kmeno1_x,U_kmeno1_y-W,k,r);
                elseif(k<=2&&r<=2)
                    [ Ax(r,k), ...
                        bx(r) ]= ...
                        AssembleVel1DXREL(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                        quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',Data.dt,Data.nu,Fx,U_kmeno1_x,U_kmeno1_y-W,k,r);
                elseif(k<=2&&r>2)
                    [ Ax(1+(r-1-2)*fedof +2: (r-2)*fedof+2 ,k), ...
                        bx(1+(r-1-2)*fedof : (r-2)*fedof ) ]= ...
                        AssembleVel1DXREL(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                        quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',Data.dt,Data.nu,Fx,U_kmeno1_x,U_kmeno1_y-W,k,r);
                elseif(k>2&&r<=2)
                    [ Ax(r , 1+(k-1-2)*fedof +2: (k-2)*fedof+2), ...
                        bx(1+(r-1)*fedof : r*fedof ) ]= ...
                        AssembleVel1DXREL(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                        quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',Data.dt,Data.nu,Fx,U_kmeno1_x,U_kmeno1_y-W,k,r);
                end
            end
        end
        
        for k = 1:mx+2
            
            DMbKX=DModalBasisX(:,k);
            
            for r = 1:my
                
                MbRY=ModalBasisY(:,r);
                DMbRY=DModalBasisY(:,r);
                if(k>2)
                    Axy(1+(r-1)*fedof : r*fedof , 1+(k-1-2)*fedof+2 : (k-2)*fedof+2)=...
                        AssembleVelXYREL(Map,weights_yhat,...
                        DMbKX,MbRY,DMbRY,...
                        FemBasisU',DFemBasisU',ne,...
                        quad_weights_mesh_x,quad.nx,Data.nu,k);
                elseif(k<=2)
                    Axy(1+(r-1)*fedof : r*fedof ,k)=...
                        AssembleVelXYREL(Map,weights_yhat,...
                        DMbKX,MbRY,DMbRY,...
                        FemBasisU',DFemBasisU',ne,...
                        quad_weights_mesh_x,quad.nx,Data.nu,k);
                end
                
            end
        end
        
        for k = 1:my
            
            MbKY=ModalBasisY(:,k);
            DMbKY=DModalBasisY(:,k);
            
            MBrobKY=MBrobY(:,k);
            
            for r = 1:my
                MBrobRY=MBrobY(:,r);
                
                MbRY=ModalBasisY(:,r);
                DMbRY=DModalBasisY(:,r);
                
                [ Ay(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof), ...
                    by(1+(r-1)*fedof : r*fedof ) ]= ...
                    FSIAssembleVel1DY(Map,weights_yhat, MbKY, MbRY,DMbKY,DMbRY,ne,...
                    quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',Data.dt,Data.nu,...
                    Fy,U_kmeno1_x,U_kmeno1_y-W,Data.K,MBrobKY,MBrobRY,Etadown_old,Etaup_old);
            end
        end
        
        for k = 1:my
            
            MbKY=ModalBasisY(:,k);
            DMbKY=DModalBasisY(:,k);
            
            for r = 1:mx+2
                
                DMbRX=DModalBasisX(:,r);
                if(r<=2)
                    Ayx(r, 1+(k-1)*fedof : k*fedof)=...
                        AssembleVelYXREL(Map,weights_yhat,...
                        MbKY,DMbKY,DMbRX,...
                        FemBasisU',DFemBasisU',ne,...
                        quad_weights_mesh_x,quad.nx,Data.nu,r);
                elseif(r>2)
                    Ayx(1+(r-1-2)*fedof+2 : (r-2)*fedof +2, 1+(k-1)*fedof : k*fedof)=...
                        AssembleVelYXREL(Map,weights_yhat,...
                        MbKY,DMbKY,DMbRX,...
                        FemBasisU',DFemBasisU',ne,...
                        quad_weights_mesh_x,quad.nx,Data.nu,r);
                end
            end
        end
        
        for k = 1:mx+2
            MbKX=ModalBasisX(:,k);
            DMbKX=DModalBasisX(:,k);
            
            for s=1:mp
                
                MbSP=ModalBasisP(:,s);
                if(k>2)
                    Bx(1+(s-1)*fedofP:s*fedofP, 1+(k-1-2)*fedof+2 : (k-2)*fedof+2) = ...
                        AssembleBXREL(Map,weights_yhat,MbKX,MbSP,DMbKX,ne,quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',FemBasisP',k);
                elseif(k<=2)
                    Bx(1+(s-1)*fedofP:s*fedofP, k) = ...
                        AssembleBXREL(Map,weights_yhat,MbKX,MbSP,DMbKX,ne,quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',FemBasisP',k);
                    
                end
            end
        end
        
        for k = 1:my
            
            DMbKY=DModalBasisY(:,k);
            
            for s=1:mp
                
                MbSP=ModalBasisP(:,s);
                
                By(1+(s-1)*fedofP:s*fedofP, 1+(k-1)*fedof : k*fedof) = ...
                    AssembleBY(Map,weights_yhat,MbSP,DMbKY,ne,quad_weights_mesh_x,quad.nx,FemBasisU',FemBasisP');
            end
        end
                
        %[Ay,by]=DirichletInflow(Ay,by,zeros(my,1),my,fedof); %uy=0 all'inflow
        %[Ay,by]=DirichletOutflow(Ay,by,zeros(my,1),my,fedof); %uy=0 all'outflow
        %BC per u in inflow
        if(Data.essential_inflow)
            [Ax,bx]=DirichletInflowREL(Ax,bx,Data.inflow_dirichlet,mx+2,fedof);
        else
            bx=NeumannInflowREL(weights_yhat,bx,Data.Pinf(t),ModalBasisX,fedof,Jin_out);
        end
        bx=NeumannOutflowREL(weights_yhat,bx,Data.Pout(t),ModalBasisX,fedof,Jin_out);
        
        
        % FINE ASSEMBLAGGIO %
        A=[Ax,Ayx,Bx';Axy,Ay,By';Bx,By,zeros(ndofP,ndofP)];

        % Soluzione sistema lineare
        sol=A\[bx;by;zeros(ndofP,1)];
        u_k_x=sol(1:ndofx+2);
        u_k_y=sol(ndofx+1+2:ndofx+ndofy+2);
        
        %Calcola il nuovo spostamento
        [etadown_kp1,etaup_kp1]=ComputeDisplacement(u_k_y,Data.dt,etadown_old,etaup_old,MBrobY);
        %Calcola tutti i nuovi Jacobiani e simili
        [Map,Y,W,Ynodes,Jin_out]=ComputeMapFromDisplacement(etadown_kp1,etaup_kp1,Xhat,Yold,Yhat,FemBasisU,DFemBasisU,Data.dt);
                                                            
        test=sqrt(sum((etadown_kp1-etadown_k).^2+(etaup_kp1-etaup_k).^2)) /sqrt(sum((etadown_kp1).^2+(etaup_kp1).^2));
        test=test+sqrt(sum((u_kmeno1_x-u_k_x).^2)+sum((u_kmeno1_y-u_k_y).^2))/sqrt(sum(u_k_x.^2)+sum(u_k_y.^2))
        etadown_k=etadown_kp1;
        etaup_k=etaup_kp1;
        u_kmeno1_x=u_k_x;
        u_kmeno1_y=u_k_y;
    end
	Yold=Y;
    display(['Ended fixed point iterations Nb = ', num2str(PF_iter), ' max= ' , num2str(PF_max)])
    display(['min displacement = ',num2str(min(min(-etadown_k),min(etaup_k))),'max displacement = ',num2str(max(max(-etadown_k),max(etaup_k)))])
    uold_x=sol(1:ndofx+2);
    uold_y=sol(ndofx+1+2:ndofx+ndofy+2);
    
    %Assegna lo spostamento etaold
    etadown_old=etadown_k;
    etaup_old=etaup_k;
    % Post Processing
    if (t_iter==Data.ntimestep)
        out.ended=true;
    end
    plot_stokesREL(sol,mesh_x,Ynodes,ModalBasisX,ModalBasisY,ModalBasisP,out,t,etadown_old,etaup_old,mx,my,mp);
end
% A=[Ax,Axy,Bx';Ayx,Ay,By';Bx,By,zeros(ndofP,ndofP)];
% figure;spy(A);
