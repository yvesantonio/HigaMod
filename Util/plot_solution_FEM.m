function [errL2,errH1,aus_figure,sol,meshx,verMesh] = plot_solution_FEM(size_mb,a_ril,b_ril,cutx,...
                                                                              hx, u,bc_up,bc_down,Coeff_forma,...
                                                                              caso,sol_exx,forma,dforma,psJac,...
                                                                              psJac2,p,k,soldx,soldy,geometricInfo)

    % SETTING FOR THE PERSISTENT CALLS OF THE FIGURE PLOT

    persistent calls current_figure

    if isempty(calls)
        calls = 0;
        current_figure = figure;
    else
        figure(current_figure);
    end	

    if calls < 10
        strcalls = strcat('000',int2str(calls));
    elseif calls < 100
        strcalls = strcat('00',int2str(calls));
    elseif calls >= 100
        strcalls = strcat('0',int2str(calls));
    end

    len = sum(psJac);
    colormap jet;

    % SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    numbNodes = ne + 1;

    % SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 64;

    % CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

    verMesh = linspace(0,1,M);

    y_elements = length(verMesh);

    % INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(numbNodes,y_elements);
    sol_x = zeros(numbNodes,y_elements);
    sol_y = zeros(numbNodes,y_elements);

    import Core.BasisHandler

    obj_newModalBasis = BasisHandler();

    obj_newModalBasis.dimModalBasis         = size_mb;
    obj_newModalBasis.evalNodesY            = verMesh;
    obj_newModalBasis.labelUpBoundCond      = bc_up{1};
    obj_newModalBasis.labelDownBoundCond    = bc_down{1};
    obj_newModalBasis.coeffForm             = Coeff_forma;

    [coeffModalBase,coeffModalBaseDer] = newModalBasis(obj_newModalBasis);

    %---------------------------------------------------------------------%
    % Note:
    % The coefficients 'coeffModalBase' and 'coeffModalBaseDer' are
    % vectors conrresponding respectively to the coefficients of the
    % modal bases and derivative of the modal basis on the points
    % assigned in the Y direction.
    %---------------------------------------------------------------------%

    % SETTING OF THE FINITE ELEMENT MESH IN THE X DIRECTION

    %---------------------------------------------------------------------%
    % NOTE:
    % Here we are considering a linear parametrization of the problem
    % and the continuity requirements achieved are given by C^(p-k).
    %---------------------------------------------------------------------%

    horMesh = linspace(0,1*len,numbNodes);
    
     in = 0; fi = 1*len;
     knot      = augknt([in fi],p+1);
     h         = hx;
     nel       = ne;
     ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
     ins       = (fi-in)*ins+in;
     [~,knot] = bspkntins(p,in:(fi-in)*1/p:fi,knot,ins);
     % xx2 = linspace(0,1*len,2*localdata_upBDomain);
    
    % SETTING OF THE FINITE ELEMENT MESH IN THE X DIRECTION
    
    [X,Y] = meshgrid(horMesh,verMesh);

    L_computed = geometricInfo.L(X);
    a_computed = geometricInfo.a(X);
    
    Y = L_computed.*Y+a_computed;
    
    % EVALUATION OF THE SOLUTION
    
    for h=1:numbNodes
        for k=1:M
            for imb=1:size_mb
                sol(h,k)=sol(h,k)+ u(h+(imb-1)*numbNodes)*coeffModalBase(k,imb);
            end
        end
    end

%     for m = 1:size_mb
%         approx((m-1)*2*nnx+1:2*m*nnx)  = bspeval(p,u((m-1)*nnx+1:m*nnx)',knot,xx2);
%         [dc,dk] = bspderiv(p,u((m-1)*nnx+1:m*nnx)',knot);
%         dapprox((m-1)*2*nnx+1:2*m*nnx) = bspeval(p-1,dc,dk,xx2);
%         for j = 2:2*nnx
%              dapprox((m-1)*2*nnx+j) = dapprox((m-1)*2*nnx+j);
%         end
%     end
%      
%     u = approx';
%     solsol = approx';
%     
%     % COMPUTATION OF THE NURBS BASE
% 
%     u  = approx';
% 
%     % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE
% 
%     ux = dapprox';
% 
%     solsol = approx';
% 
%     for h = 1:2*nnx
% 
%         for k = 1:M
% 
%             for imb = 1:size_mb
% 
%                 % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
%                 % EVALUATED IN THE POINTS OF THE DOMAIN MESH
% 
%                 sol(h,k)   = sol(h,k)   + u(h+(imb-1)*2*nnx)*coeffModalBase(k,imb);
% 
%                 % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
%                 % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
%                 % DOMAIN MESH
% 
%                 sol_x(h,k) = sol_x(h,k) + ux(h+(imb-1)*2*nnx)*coeffModalBase(k,imb);
%                 sol_y(h,k) = sol_y(h,k) + u(h+(imb-1)*2*nnx)*coeffModalBaseDer(k,imb);
% 
%             end
% 
%         sol(h,k) = sol(h,k)+ a_ril(1) * verMesh(k) + b_ril(1);
% 
%         end
% 
%     end

    %-----------------------------------------------------------------%
    % NOTE:
    % Up to this point we have computed the matrix containing the
    % approximated solution of the differential problem. However, in
    % order to plot correctly the vector field we must change the point
    % defining the domain as follows. We will create the new matrix X
    % and Y, with the same size, to correctly represent our domain.
    %-----------------------------------------------------------------%

    leftLim     = 0.5;
    rightLim    = 0.5;
    
    X = zeros(M,numbNodes);
    Y = zeros(M,numbNodes);
    k = linspace(cutx(1),cutx(end),numbNodes);
    inv = 1;

    for kk = 1:numbNodes

        pp  = k(kk);

        if kk > 1

            pp1 = k(kk-1);

            if ((-1/dforma(pp))*(-1/dforma(pp1))<0)
                inv = inv + 1;
            end

        end

        xx = linspace(pp+(-1)^(inv)*leftLim*sin(-pi/2+atan(-1/dforma(pp))),pp-(-1)^(inv)*rightLim*sin(-pi/2+atan(-1/dforma(pp))),M);
        yy = linspace(forma(pp)-(-1)^(inv)*leftLim*cos(-pi/2+atan(-1/dforma(pp))),forma(pp)+(-1)^(inv)*rightLim*cos(-pi/2+atan(-1/dforma(pp))),M);

        X(:,kk) = xx;
        Y(:,kk) = yy;

    end

    minX = min(min(min(X),0));
    maxX = max(max(max(X),0));
    minY = min(min(min(Y),0));
    maxY = max(max(max(Y),0));
    
    %---------------------------------------------------------------------%
    %       CONTOUR PLOT OF THE APPROXIMATED SOLUTION IN THE DOMAIN
    %---------------------------------------------------------------------%
    figure;
            
    [~,~] = contourf(X,Y,sol',20);    
    colormap(jet);
    cmin = min(min(sol'));
    cmax = max(max(sol'));
    caxis([cmin cmax])
    colorbar();
    %axis([minX maxX minY maxY]);
    %axis equal
    set(gca, 'FontSize', 14)

    %--------------------------------------------------%
    % CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
    %--------------------------------------------------%

    %---------------------------------------------------------------------%
    % Note:
    % The following lines of code represent the plot of the exact solution
    % of the ADR problem considering only the case number 1. More testing
    % files will be added later on.
    %---------------------------------------------------------------------%

    set(gca,'FontUnits','points','FontSize',11);
    colorbar;

    calls = calls + 1;
    aus_figure = current_figure;

    size_fb     = numbNodes;
    meshx       = linspace(0,cutx(end),size_fb);
    verMesh  = linspace(-0.5,0.5,M);

    %---------------------------------------------------------------------%
    %                   EVALUATION OF THE SOLUTION ERROR
    %---------------------------------------------------------------------%
    
    sol_ex      = zeros(M,numbNodes);
    sol_exdx    = zeros(M,numbNodes);
    sol_exdy    = zeros(M,numbNodes);

%     errL2 = 0;
%     errH1 = 0;
%     aus_figure = 0;
%     solsol = 0;
%     meshx = 0;
%     evalNodesY = 0;
    
    % CREATION OF THE MESH TO EVALUATE THE EXACT SOLUTION OF THE PROBLEM 
    
    exactMeshX = zeros(1,numbNodes);
    exactMeshX(1) = 0;

    for ii = 2 : numbNodes
        exactMeshX(ii) = exactMeshX(ii - 1) + psJac(ii - 1);
    end
    
    exactMeshY = linspace(-0.5,0.5,M);

    for ii = 1:M
        
        for j = 1:numbNodes
            
            % SETTING THE POINTS TO EVALUATE THE EXACT SOLUTION
            
            xv          = exactMeshX(j); 
            yv          = exactMeshY(ii);
            
            % EVALUATION OF THE EXACT SOLUTION
            
            sol_ex(ii,j)   = sol_exx(xv,yv);
            sol_exdx(ii,j) = soldx(xv,yv);
            sol_exdy(ii,j) = soldy(xv,yv);
        end
        
    end
    
    % EVALUATION OF THE ERROR OF THE SOLUTION AND ITS DERIVATIVE IN THE X
    % AND Y DIRECTIONS
    
    err1 = abs(sol_ex'   - sol  );
    err2 = abs(sol_exdx' - sol_x);
    err3 = abs(sol_exdy' - sol_y);
    
    % COMPUTATION OF THE L2 NORM OF THE ERROR
    
    errL2 = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2)));
    
    % COMPUTATION OF THE H1 NORM OF THE ERROR
    
    errH1 = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2))  +...
                 sum(sum(err2.^2))/sum(sum(sol_exdx.^2))+...
                 sum(sum(err3.^2))/sum(sum(sol_exdy.^2)));

end
