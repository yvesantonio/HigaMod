function [errL2,errH1,solsol,meshx,evalNodesY] = plot_solution_IGA(size_mb,a_ril,b_ril,cutx,...
                                                                              hx, u,bc_up,bc_down,Coeff_forma,...
                                                                              caso,sol_exx,forma,dforma,psJac,...
                                                                              p,k,soldx,soldy,cuty)

%     [HbFunc,HmFunc] = observStates(size_mb,cutx,hx,bc_up,bc_down,Coeff_forma,p,k,cuty);

    % SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne * k + p + 1 - k;           % Number of Control Point on the X Direction
                                        % for the IGA Mesh

    % SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 64;

    % CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

    evalNodesY = linspace(cuty(1),cuty(2),M);

    y_elements = length(evalNodesY);

    % INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(nnx,y_elements);
    sol_x = zeros(nnx,y_elements);
    sol_y = zeros(nnx,y_elements);

    import Core.BasisHandler

    obj_newModalBasis = BasisHandler();

    obj_newModalBasis.dimModalBasis         = size_mb;
    obj_newModalBasis.evalNodesY            = evalNodesY;
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

    % SETTING OF THE ISOGEOMETRIC MESH IN THE X DIRECTION

    %---------------------------------------------------------------------%
    % NOTE:
    % Here we are considering a linear parametrization of the problem
    % and the continuity requirements achieved are given by C^(p-k).
    %---------------------------------------------------------------------% 

    knot      = augknt([cutx(1) cutx(2)],p+1);
    h         = hx;
    internKnot = linspace(0,1,ne+1)';
    internKnot = internKnot(2:end-1);
    ins       = sort(reshape(internKnot*ones(1,k),1,[]));%k*(nel-1)));
    ins       = (cutx(2)-cutx(1))*ins + cutx(1);

    [~,knot] = bspkntins(p,cutx(1):(cutx(2)-cutx(1))*1/p:cutx(2),knot,ins);

    %---------------------------------------------------------------------%
    % Note:
    % The vector 'knot' contains the nodes used in the mesh.
    %---------------------------------------------------------------------%
    
    xx3 = linspace(cutx(1),cutx(2),nnx)';

    %---------------------------------------------------------------------%
    % Note:
    % The function 'bspeval' used in the following loops is taken from
    % the NURBS toolbox and its primary objective is to evaluate a 
    % univariate B-Spline function in the desired nodes.
    %
    % The parameter involved in this functions are (in order):
    %
    % (1ST INPUT) Degree of the B-Spline;
    % (2ND INPUT) Control Points, matrix of size (dim,nc), where 'dim' 
    %             is the degree of the B-Spline function and 'nc' is
    %             the number of control points.
    % (3RD INPUT) Knot sequence, row vector of size nk, where 'nk' is
    %             the number of knot points.
    % (4TH INPUT) Parametric evaluation points, row vector of size nu,
    %             where 'nu' is the number of points used to represent
    %             the domain.
    %
    % (OUTPUT) Evaluated points, matrix of size (dim,nu).
    %---------------------------------------------------------------------%

    Hb = [];
    
    solsol = u;
    
    switch caso
            case {1,2,3,4,5,6,7,8,9,10}
            for m = 1:size_mb
                
                approx((m-1)*nnx+1:m*nnx)  = bspeval(p,u((m-1)*nnx+1:m*nnx)',knot,xx3);
                
                % CREATION OF THE BASIS MATRIX
                
                Hb_temp{m} = basisMatrix(p,u((m-1)*nnx+1:m*nnx)',knot,xx3);
                Hb = blkdiag(Hb,Hb_temp{m});
                
                [dc,dk] = bspderiv(p,u((m-1)*nnx+1:m*nnx)',knot);
                dapprox((m-1)*nnx+1:m*nnx) = bspeval(p-1,dc,dk,xx3);

                for j = 2:nnx
                     dapprox((m-1)*nnx+j) = dapprox((m-1)*nnx+j)/sqrt(1+(dforma(xx3(j)))^2);
                end

            end

            % COMPUTATION OF THE NURBS BASE

            u  = approx';

            % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE

            ux = dapprox';
    end
    
    Hm = zeros(nnx * M, nnx * size_mb);
    solVect = zeros(nnx * M,1);
    solMat = zeros(M,nnx);

    for h = 1:nnx

        for k = 1:M

            for imb = 1:size_mb
                
                % CREATION OF THE MODAL BASIS MATRIX
                
                Hm(h + (k - 1) * nnx , h + (imb - 1) * nnx) = coeffModalBase(k,imb);

                % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
                % EVALUATED IN THE POINTS OF THE DOMAIN MESH

                sol(h,k)   = sol(h,k)   + u(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                solVect(h + (k-1)*nnx) = solVect(h + (k-1)*nnx) + u(h+(imb-1)*nnx)*coeffModalBase(k,imb);

                % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
                % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
                % DOMAIN MESH

                sol_x(h,k) = sol_x(h,k) + ux(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                sol_y(h,k) = sol_y(h,k) + u(h+(imb-1)*nnx)*coeffModalBaseDer(k,imb);

            end

        sol(h,k) = sol(h,k)+ a_ril(1) * evalNodesY(k) + b_ril(1);

        end

    end
    
    finalSol = Hm * Hb * solsol;
    for ii = 1:M
        solMat(ii,:) = finalSol((ii-1)*nnx + 1 : ii * nnx);
    end
    
%     disp('Matrix Errors')
%     disp(norm(Hm - HmFunc,inf));
%     disp(norm(Hb - HbFunc,inf));
%     
%     disp('Solution error')
%     disp(norm(solMat - sol',inf));

    %-----------------------------------------------------------------%
    % NOTE:
    % Up to this point we have computed the matrix containing the
    % approximated solution of the differential problem. However, in
    % order to plot correctly the vector field we must change the point
    % defining the domain as follows. We will create the new matrix X
    % and Y, with the same size, to correctly represent our domain.
    %-----------------------------------------------------------------%

    downLim  = cuty(1);
    upLim    = cuty(2);
    
    X = zeros(M,nnx);
    Y = zeros(M,nnx);
    k = linspace(cutx(1),cutx(end),nnx);
    inv = 1;

    for kk = 1:nnx

        pp  = k(kk);

        if kk > 1

            pp1 = k(kk-1);

            if ((-1/dforma(pp))*(-1/dforma(pp1))<0)
                inv = inv + 1;
            end

        end
        
        % ATTENTION!!!!
        %-----------------------------------------------------------------%
        % The original code in the next two lines work only when the point
        % in the Y direction are symmetric with respect to the X axis. I
        % perform some modifications to change this, but it still does not
        % work for the most general case (channel thickness varying).
        %-----------------------------------------------------------------%  

        xx = linspace(pp+(-1)^(inv)*((cuty(2) - cuty(1))/2)*sin(-pi/2+atan(-1/dforma(pp))),...
                      pp-(-1)^(inv)*((cuty(2) - cuty(1))/2)*sin(-pi/2+atan(-1/dforma(pp))),M);

        yy = linspace(((cuty(2) + cuty(1))/2) + forma(pp)-(-1)^(inv)*((cuty(2) - cuty(1))/2)*cos(-pi/2+atan(-1/dforma(pp))),...
                       ((cuty(2) + cuty(1))/2) + forma(pp)+(-1)^(inv)*((cuty(2) - cuty(1))/2)*cos(-pi/2+atan(-1/dforma(pp))),M);

        X(:,kk) = xx;
        Y(:,kk) = yy;

    end
    
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    
%     [p,e,t,ffSol] = my_FFimportfilemesh('ffmesh.msh','ffSolution.sol');
%     SOL = griddata(p(1,:),p(2,:),ffSol,X,Y,'linear');
%     
%     disp(SOL);
%     
%     figure;
%     pdeplot(p,e,t,'xydata',ffSol,'mesh','on','colormap','jet');
%     
%     figure;
%     mesh(X,Y,SOL);
    
    
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%

    
    
%     figure;
            
%     [~,~] = contourf(X,Y,SOL,20);    
%     colormap(parula);
%     cmin = min(min(solMat));%-0.05;%min(min(sol'));
%     cmax = max(max(solMat)); %0.2;% max(max(sol'));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 14)

    %--------------------------------------------------%
    % CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
    %--------------------------------------------------%

    %---------------------------------------------------------------------%
    % Note:
    % The following lines of code represent the plot of the exact solution
    % of the ADR problem considering only the case number 1. More testing
    % files will be added later on.
    %---------------------------------------------------------------------%

    size_fb     = nnx;
    meshx       = linspace(cutx(1),cutx(2),size_fb);
    evalNodesY  = linspace(cuty(1),cuty(2),M);

    %---------------------------------------------------------------------%
    %                   EVALUATION OF THE SOLUTION ERROR
    %---------------------------------------------------------------------%
    
    sol_ex      = zeros(M,nnx);
    sol_exdx    = zeros(M,nnx);
    sol_exdy    = zeros(M,nnx);
    
    % CREATION OF THE MESH TO EVALUATE THE EXACT SOLUTION OF THE PROBLEM 
    
%     exactMeshX = zeros(1,nnx);
%     exactMeshX(1) = cutx(1);
% 
%     for ii = 2 : nnx
% %         exactMeshX(ii) = exactMeshX(ii - 1) + psJac(ii - 1);
%         exactMeshX(ii) = exactMeshX(ii - 1) + psJac(ii - 1);
%     end
    
    exactMeshX = linspace(cutx(1),cutx(2),nnx);
    exactMeshY = linspace(cuty(1),cuty(2),M);

    for ii = 1:M
        
        for j = 1:nnx
            
            % SETTING THE POINTS TO EVALUATE THE EXACT SOLUTION
            
            xv          = exactMeshX(j); 
            yv          = exactMeshY(ii);
            
            % EVALUATION OF THE EXACT SOLUTION
            
            sol_ex(ii,j)   = sol_exx(xv,yv);
            sol_exdx(ii,j) = soldx(xv,yv);
            sol_exdy(ii,j) = soldy(xv,yv);
        end
        
    end
    
    %---------------------------------------------------------------------%
    %          CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
    %---------------------------------------------------------------------%
%     figure;
%             
%     [~,~] = contourf(exactMeshX,exactMeshY,sol_ex,20);    
%     colormap(jet);
%     cmin = min(min(sol_ex));
%     cmax = max(max(sol_ex));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 14)
    
    % EVALUATION OF THE ERROR OF THE SOLUTION AND ITS DERIVATIVE IN THE X
    % AND Y DIRECTIONS
    
    err1 = abs(sol_ex'   - sol  );
    err2 = abs(sol_exdx' - sol_x);
    err3 = abs(sol_exdy' - sol_y);
    
    % COMPUTATION OF THE L2 NORM OF THE ERROR
    
    errL2 = norm(err1,2)/norm(sol_ex,2);
    
    % COMPUTATION OF THE H1 NORM OF THE ERROR

    errH1 = norm(err1,2)/norm(sol_ex,2) + norm(err2,2)/norm(sol_exdx,2) + norm(err3,2)/norm(sol_exdy,2);
    
    minX = min(min(min(X),0));
    maxX = max(max(max(X),0));
    minY = min(min(min(Y),0));
    maxY = max(max(max(Y),0));
    
    %---------------------------------------------------------------------%
    %       CONTOUR PLOT OF THE APPROXIMATED SOLUTION IN THE DOMAIN
    %---------------------------------------------------------------------%
    
    figure;
            
    [~,~] = contourf(X,Y,sol_ex,20);    
    colormap(parula);
    cmin = min(min(sol'));%-0.05;%min(min(sol'));
    cmax = max(max(sol')); %0.2;% max(max(sol'));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 14)
    
    figure;
            
    [~,~] = contourf(X,Y,solMat,20);    
    colormap(parula);
    cmin = min(min(solMat));%-0.05;%min(min(sol'));
    cmax = max(max(solMat)); %0.2;% max(max(sol'));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 14)

end
