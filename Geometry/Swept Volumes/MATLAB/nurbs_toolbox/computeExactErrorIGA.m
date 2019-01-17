function [errL2,errH1] = computeExactErrorIGA(size_mb,a_ril,b_ril,cutx,...
                                          hx, u,bc_up,bc_down,Coeff_forma,...
                                          caso,sol_exx,forma,dforma,psJac,...
                                          p,horNodes,soldx,soldy,cuty)

%     [HbFunc,HmFunc] = observStates(size_mb,cutx,hx,bc_up,bc_down,Coeff_forma,p,k,cuty);

    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

    % SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    numbHorNodes = ne * horNodes + p + 1 - horNodes;

    % CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS
    
    % Compute the Gauss-Legendre Nodes
    
    numbQuadNodes = 100;
    numbEvalNodes = (numbHorNodes-1) * numbQuadNodes;
    M = numbQuadNodes;
    
    % QUADRATURE RULE FOR THE HORIZONTAL DIRECTION

    objGL = IntegrateHandler();
    objGL.numbQuadNodes = numbQuadNodes;
    [~, refGLNodes, refGLWeights] = gaussLegendre(objGL);
    
    % QUADRATURE RULE FOR THE VERTICAL DIRECTION

    evalNodesY = linspace(cuty(1),cuty(2),M);
    verWeights = weiCompSimpson(numbQuadNodes,cuty(1),cuty(2));

    y_elements = length(evalNodesY);

    % INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(numbEvalNodes,y_elements);
    sol_x = zeros(numbEvalNodes,y_elements);
    sol_y = zeros(numbEvalNodes,y_elements);

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
    internKnot = linspace(0,1,ne + 1)';
    internKnot = internKnot(2:end-1);
    ins       = sort(reshape(internKnot*ones(1,horNodes),1,[]));
    ins       = (cutx(2)-cutx(1))*ins + cutx(1);

    [~,knot] = bspkntins(p,cutx(1):(cutx(2)-cutx(1))*1/p:cutx(2),knot,ins);

    %---------------------------------------------------------------------%
    % Note:
    % The vector 'knot' contains the nodes used in the mesh.
    %---------------------------------------------------------------------%
    
    hNodes = linspace(cutx(1),cutx(2),numbHorNodes);
    
    for i = 1:(numbHorNodes-1)
        obj_quadratureRule = IntegrateHandler();
                
        obj_quadratureRule.leftBoundInterval = hNodes(i);
        obj_quadratureRule.rightBoundInterval = hNodes(i+1);
        obj_quadratureRule.inputNodes = refGLNodes;
        obj_quadratureRule.inputWeights = refGLWeights;

        [horNodes((i-1)*numbQuadNodes+1 : i*numbQuadNodes), ...
         horWeights((i-1)*numbQuadNodes+1 : i*numbQuadNodes)] = ...
         quadratureRule(obj_quadratureRule);
    end

%     horNodes = linspace(cutx(1),cutx(2),numbEvalNodes)';
    
    % DEBUG
    figure
    plot(horNodes);
    hold on;
    plot(horWeights);

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

    % Hb = [];
    
    % solsol = u;
    
    switch caso
            case {1,2,3,4,5,6,7,8,9,10}
            for m = 1:size_mb
                
                approx((m-1)*numbEvalNodes+1:m*numbEvalNodes)  = bspeval(p,u((m-1)*numbHorNodes+1:m*numbHorNodes)',knot,horNodes);
                
                % CREATION OF THE BASIS MATRIX
                
                % Hb_temp{m} = basisMatrix(p,u((m-1)*numbHorNodes+1:m*numbHorNodes)',knot,horNodes);
                % Hb = blkdiag(Hb,Hb_temp{m});
                
                [dc,dk] = bspderiv(p,u((m-1)*numbHorNodes+1:m*numbHorNodes)',knot);
                dapprox((m-1)*numbEvalNodes+1:m*numbEvalNodes) = bspeval(p-1,dc,dk,horNodes);

                for j = 2:numbEvalNodes
                     dapprox((m-1)*numbEvalNodes+j) = dapprox((m-1)*numbEvalNodes+j)/sqrt(1+(dforma(horNodes(j)))^2);
                end

            end

            % COMPUTATION OF THE NURBS BASE

            u  = approx';

            % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE

            ux = dapprox';
    end
    
    % Hm = zeros(numbEvalNodes * M, numbEvalNodes * size_mb);
    solVect = zeros(numbEvalNodes * M,1);

    for horPos = 1:numbEvalNodes

        for verPos = 1:M

            for imb = 1:size_mb
                
                % CREATION OF THE MODAL BASIS MATRIX
                
                % Hm(horPos + (verPos - 1) * numbEvalNodes , horPos + (imb - 1) * numbEvalNodes) = coeffModalBase(verPos,imb);

                % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
                % EVALUATED IN THE POINTS OF THE DOMAIN MESH

                sol(horPos,verPos)   = sol(horPos,verPos)   + u(horPos+(imb-1)*numbEvalNodes)*coeffModalBase(verPos,imb);
                solVect(horPos + (verPos-1)*numbEvalNodes) = solVect(horPos + (verPos-1)*numbEvalNodes) + u(horPos+(imb-1)*numbEvalNodes)*coeffModalBase(verPos,imb);

                % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
                % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
                % DOMAIN MESH

                sol_x(horPos,verPos) = sol_x(horPos,verPos) + ux(horPos+(imb-1)*numbEvalNodes)*coeffModalBase(verPos,imb);
                sol_y(horPos,verPos) = sol_y(horPos,verPos) + u(horPos+(imb-1)*numbEvalNodes)*coeffModalBaseDer(verPos,imb);

            end

        sol(horPos,verPos) = sol(horPos,verPos)+ a_ril(1) * evalNodesY(verPos) + b_ril(1);

        end

    end
    
    solMat = sol;
    
    %-----------------------------------------------------------------%
    % NOTE:
    % Up to this point we have computed the matrix containing the
    % approximated solution of the differential problem. However, in
    % order to plot correctly the vector field we must change the point
    % defining the domain as follows. We will create the new matrix X
    % and Y, with the same size, to correctly represent our domain.
    %-----------------------------------------------------------------%
    
    X = zeros(M,numbEvalNodes);
    Y = zeros(M,numbEvalNodes);
    horNodes = linspace(cutx(1),cutx(end),numbEvalNodes);
    inv = 1;

    for kk = 1:numbEvalNodes

        pp  = horNodes(kk);

        if kk > 1

            pp1 = horNodes(kk-1);

            if ((-1/dforma(pp))*(-1/dforma(pp1))<0)
                inv = inv + 1;
            end

        end

        xx = linspace(pp+(-1)^(inv)*((cuty(2) - cuty(1))/2)*sin(-pi/2+atan(-1/dforma(pp))),...
                      pp-(-1)^(inv)*((cuty(2) - cuty(1))/2)*sin(-pi/2+atan(-1/dforma(pp))),M);

        yy = linspace(((cuty(2) + cuty(1))/2) + forma(pp)-(-1)^(inv)*((cuty(2) - cuty(1))/2)*cos(-pi/2+atan(-1/dforma(pp))),...
                      ((cuty(2) + cuty(1))/2) + forma(pp)+(-1)^(inv)*((cuty(2) - cuty(1))/2)*cos(-pi/2+atan(-1/dforma(pp))),M);

        X(:,kk) = xx;
        Y(:,kk) = yy;

    end
    
    % CREATE WEIGHT MATRIX
    %---------------------------------------------------------------------%
    % Note: We use the one-dimensional gauss-legendre quadrature nodes to
    % create the weight matrix to be used in the computation of the
    % integral representing the error.
    %---------------------------------------------------------------------%
    
    % DEBUG
    disp(size(verWeights))
    figure
    plot(verWeights);
    
    [n,m] = size(X);
    
    for ii = 1:n
        for jj = 1:m
            W(ii,jj) = verWeights(ii) * horWeights(jj);
        end
    end
    
    figure
    mesh(X,Y,W);
    
    %---------------------------------------------------------------------%
    % Compute the approximated error
    %---------------------------------------------------------------------%
    
    exactSol  = sol_exx(X,Y);
    exactSolx = soldx(X,Y);
    exactSoly = soldy(X,Y);
    
    absError  = abs(exactSol - solMat');
    absErrorx = abs(exactSolx - sol_x');
    absErrory = abs(exactSoly - sol_y');
    
    errL2 = 0;
    errH1 = 0;
    
    for ii = 1:n
        for jj = 1:m
            errL2 = errL2 + absError(ii,jj) * W(ii,jj);
            errH1 = errH1 + absErrorx(ii,jj) * W(ii,jj);
            errH1 = errH1 + absErrory(ii,jj) * W(ii,jj);
        end
    end
    
    errH1 = errH1 + errL2;
    
    %---------------------------------------------------------------------%
    % Contour plots
    %---------------------------------------------------------------------%
    
    figure;
            
    [~,~] = contourf(X,Y,solMat',20);    
    colormap(parula);
    cmin = min(min(sol'));%-0.05;%min(min(sol'));
    cmax = max(max(sol')); %0.2;% max(max(sol'));
    caxis([cmin cmax])
    colorbar();
    axis equal
    set(gca, 'FontSize', 14);
    
    figure;
            
    [~,~] = contourf(X,Y,exactSol,20);    
    colormap(parula);
    cmin = min(min(sol'));%-0.05;%min(min(sol'));
    cmax = max(max(sol')); %0.2;% max(max(sol'));
    caxis([cmin cmax])
    colorbar();
    axis equal
    set(gca, 'FontSize', 14);
    
    figure;
            
    [~,~] = contourf(X,Y,abs(exactSol - solMat'),20);    
    colormap(parula);
    cmin = min(min(abs(exactSol - solMat')));
    cmax = max(max(abs(exactSol - solMat')));
    caxis([cmin cmax])
    colorbar();
    axis equal
    set(gca, 'FontSize', 14);
end
