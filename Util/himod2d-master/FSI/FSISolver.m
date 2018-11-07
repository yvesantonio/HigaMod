function sol=FSISolver(mesh,quad,Data,out)

% structuring the mesh
ne=mesh.nb_elem;
nx=ne+1;
mesh_x=linspace(mesh.I(1),mesh.I(2),nx);
hx=mesh_x(2)-mesh_x(1);
node_xp2=mesh.I(1):hx/2:mesh.I(2);

%Setting up a quadrature rule for the y direction.
[~, node_y, weights_y] = gausslegendre( quad.ny );
% I adjust the nodes and weights to work on (-1,1)
node_yhat=node_y*2-1;
weights_yhat=weights_y*2;

% Setting up a quadrature rule for the x direction.
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
fedof=nx+ne;
fedofP=nx;

mx=Data.m-2;
my=Data.m;
mp=Data.m-2;

ndofx=mx*fedof;
ndofy=my*fedof;
ndofP=mp*fedofP;

%Calcola le basi FEM
[ FemBasisU , DFemBasisU ]=fembasis(2,mesh_x(1:2),quad_nodes_mesh_x(1:quad.nx));
FemBasisP = fembasis(1,mesh_x(1:2),quad_nodes_mesh_x(1:quad.nx));

%Calcola le basi Modali
[ModalBasisX,DModalBasisX]=modalbasis_legendre(mx, node_yhat, 'dir', 'dir');
[ModalBasisY,DModalBasisY]=modalbasis_legendre(my, node_yhat, 'rob', 'rob');
[MBrobY,~]=modalbasis_legendre(my, [1,-1], 'rob', 'rob');
ModalBasisP=modalbasis_legendre(mp, node_yhat, 'rob', 'rob');

%Inizializzazioni
t=0;
PF_max=50;

%Velocity
uold_x=zeros(ndofx,1);
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
[~,Ynodes]=meshgrid(node_xp2,node_yhat);
W=zeros(ne*quad.nx,quad.ny);

% Ciclo in tempo
for t_iter=1:Data.ntimestep
    
    test=1;
    PF_iter=0;
    t = t + Data.dt;
    display([ ' t = ', num2str(t),' ',num2str(t_iter),...
        '/', num2str(Data.ntimestep),...
        ])
    
    %Valuta la velocità old
    Uoldx=EvalHimod(Xhat,uold_x,ModalBasisX,FemBasisU,fedof,quad.nx);
    Uoldy=EvalHimod(Xhat,uold_y,ModalBasisY,FemBasisU,fedof,quad.nx);
    
    %Valuta lo spostamento della struttura
    [Etaup_old,Etadown_old]=EvalDisplacement(etaup_old,etadown_old,FemBasisU,Xhat(1,:),quad.nx);
    
    % Iterazioni di punto fisso
    while(test>Data.toll && PF_iter<PF_max)
        
        PF_iter=PF_iter+1;
        
        % Assemblaggio %
        Ax = sparse( ndofx, ndofx );
        Ay = sparse( ndofy, ndofy );
        
        Axy = sparse( ndofy, ndofx );
        Ayx = sparse( ndofx, ndofy );
        
        Bx = sparse( ndofP,ndofx );
        By = sparse( ndofP,ndofy );
        
        bx = zeros( ndofx ,1);
        by = zeros( ndofy ,1);
        
        %valuta la forzante
        Fx=Data.forcex(Xhat,Y,t)' + Uoldx/Data.dt;
        Fy=Data.forcey(Xhat,Y,t)' + Uoldy/Data.dt;
        
        %valuta il termine non lineare per il punto fisso
        U_kmeno1_x=EvalHimod(Xhat,u_kmeno1_x,ModalBasisX,FemBasisU,fedof,quad.nx);
        U_kmeno1_y=EvalHimod(Xhat,u_kmeno1_y,ModalBasisY,FemBasisU,fedof,quad.nx);
        
        %         display('Max ux')
        %         display(max(max(abs(U_kmeno1_x))))
        %         display('Max uy')
        %         display(max(max(abs(U_kmeno1_y))))
        
        %Assemblaggio
        
        %Ax
        for k = 1:mx
            MbKX=ModalBasisX(:,k);
            DMbKX=DModalBasisX(:,k);
            for r = 1:mx
                MbRX=ModalBasisX(:,r);
                DMbRX=DModalBasisX(:,r);
                
                [ Ax(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof), ...
                    bx(1+(r-1)*fedof : r*fedof ) ]= ...
                    AssembleVel1DX(Map,weights_yhat, MbKX, MbRX,DMbKX,DMbRX,ne,...
                    quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',1/Data.dt,...
                    Data.nu,Fx,U_kmeno1_x,U_kmeno1_y-W);
            end
        end
        
        %Axy
        for k = 1:mx
            DMbKX=DModalBasisX(:,k);
            for r = 1:my
                MbRY=ModalBasisY(:,r);
                DMbRY=DModalBasisY(:,r);
                
                Axy(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof)=...
                    AssembleVelXY(Map,weights_yhat,...
                    DMbKX,MbRY,DMbRY,...
                    FemBasisU',DFemBasisU',ne,...
                    quad_weights_mesh_x,quad.nx,Data.nu);
            end
        end
        %Ay
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
        
        %Ayx
        for k = 1:my
            MbKY=ModalBasisY(:,k);
            DMbKY=DModalBasisY(:,k);
            
            for r = 1:mx
                DMbRX=DModalBasisX(:,r);
                
                Ayx(1+(r-1)*fedof : r*fedof , 1+(k-1)*fedof : k*fedof)=...
                    AssembleVelYX(Map,weights_yhat,...
                    MbKY,DMbKY,DMbRX,...
                    FemBasisU',DFemBasisU',ne,...
                    quad_weights_mesh_x,quad.nx,Data.nu);
            end
        end
        
        %Bx
        for k = 1:mx
            MbKX=ModalBasisX(:,k);
            DMbKX=DModalBasisX(:,k);
            
            for s=1:mp
                
                MbSP=ModalBasisP(:,s);
                
                Bx(1+(s-1)*fedofP:s*fedofP, 1+(k-1)*fedof : k*fedof) = ...
                    AssembleBX(Map,weights_yhat,MbKX,MbSP,DMbKX,ne,quad_weights_mesh_x,quad.nx,FemBasisU',DFemBasisU',FemBasisP');
            end
        end
        %By
        for k = 1:my
            
            DMbKY=DModalBasisY(:,k);
            
            for s=1:mp
                
                MbSP=ModalBasisP(:,s);
                
                By(1+(s-1)*fedofP:s*fedofP, 1+(k-1)*fedof : k*fedof) = ...
                    AssembleBY(Map,weights_yhat,MbSP,DMbKY,ne,quad_weights_mesh_x,quad.nx,FemBasisU',FemBasisP');
            end
        end
        
        [Ax,Ay,bx,by]=BChandler(Ax,Ay,bx,by,Data.BC,mx,my,fedof,weights_yhat,ModalBasisX,Jin_out',t);
        % %
        
        
        % FINE ASSEMBLAGGIO %
        A=[Ax,Ayx,Bx';Axy,Ay,By';Bx,By,zeros(ndofP,ndofP)];
        
        % Soluzione sistema lineare
        sol=A\[bx;by;zeros(ndofP,1)];
        u_k_x=sol(1:ndofx);
        u_k_y=sol(ndofx+1:ndofx+ndofy);
        
        %Calcola il nuovo spostamento
        [etadown_kp1,etaup_kp1]=ComputeDisplacement(u_k_y,Data.dt,etadown_old,etaup_old,MBrobY);
        
        %         plot_stokes(sol,mesh_x,Ynodes,ModalBasisX,ModalBasisY,ModalBasisP,out,t,etadown_k,etaup_k,mx,my,mp);
        %         pause
        
        %Calcola tutti i nuovi Jacobiani e simili
        [Map,Y,W,Ynodes,Jin_out]=ComputeMapFromDisplacement(etadown_kp1,etaup_kp1,Xhat,Yold,Yhat,FemBasisU,DFemBasisU,Data.dt);
        
        test=sqrt(sum((etadown_kp1-etadown_k).^2+(etaup_kp1-etaup_k).^2)) /sqrt(sum((etadown_kp1).^2+(etaup_kp1).^2));
        test=test+sqrt(sum((u_kmeno1_x-u_k_x).^2)+sum((u_kmeno1_y-u_k_y).^2))/sqrt(sum(u_k_x.^2)+sum(u_k_y.^2));
        display(test)
        
        etadown_k=etadown_kp1;
        etaup_k=etaup_kp1;
        u_kmeno1_x=u_k_x;
        u_kmeno1_y=u_k_y;
        %         display('max W')
        %         display(max(max(abs(W))))
        
    end
    Yold=Y;
    display(['Fixed point iterations Nb = ', num2str(PF_iter)])
    if(PF_iter==PF_max )
        warning(['Maximum number of fixed point iterations reached test quantity = ',num2str(test)]);
    end
    display(['min displacement = ',num2str(min(min(-etadown_k),min(etaup_k))),' max displacement = ',num2str(max(max(-etadown_k),max(etaup_k)))])
    uold_x=sol(1:ndofx);
    uold_y=sol(ndofx+1:ndofx+ndofy);
    
    %Assegna lo spostamento etaold
    etadown_old=etadown_k;
    etaup_old=etaup_k;
    % Post Processing
    if (t_iter==Data.ntimestep)
        out.ended=true;
    end
    plot_stokes(sol,mesh_x,Ynodes,ModalBasisX,ModalBasisY,ModalBasisP,out,t,etadown_old,etaup_old,mx,my,mp);
end