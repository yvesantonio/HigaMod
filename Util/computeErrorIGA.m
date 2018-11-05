function [errL2,errH1] = computeErrorIGA(size_mb,a_ril,b_ril,cutx,...
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

    nnx = ne * horNodes + p + 1 - horNodes;           % Number of Control Point on the X Direction
                                        % for the IGA Mesh

    % SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 100;

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
    ins       = sort(reshape(internKnot*ones(1,horNodes),1,[]));%k*(nel-1)));
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

        for horNodes = 1:M

            for imb = 1:size_mb
                
                % CREATION OF THE MODAL BASIS MATRIX
                
                Hm(h + (horNodes - 1) * nnx , h + (imb - 1) * nnx) = coeffModalBase(horNodes,imb);

                % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
                % EVALUATED IN THE POINTS OF THE DOMAIN MESH

                sol(h,horNodes)   = sol(h,horNodes)   + u(h+(imb-1)*nnx)*coeffModalBase(horNodes,imb);
                solVect(h + (horNodes-1)*nnx) = solVect(h + (horNodes-1)*nnx) + u(h+(imb-1)*nnx)*coeffModalBase(horNodes,imb);

                % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
                % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
                % DOMAIN MESH

                sol_x(h,horNodes) = sol_x(h,horNodes) + ux(h+(imb-1)*nnx)*coeffModalBase(horNodes,imb);
                sol_y(h,horNodes) = sol_y(h,horNodes) + u(h+(imb-1)*nnx)*coeffModalBaseDer(horNodes,imb);

            end

        sol(h,horNodes) = sol(h,horNodes)+ a_ril(1) * evalNodesY(horNodes) + b_ril(1);

        end

    end
    
    finalSol = Hm * Hb * solsol;
    for ii = 1:M
        solMat(ii,:) = finalSol((ii-1)*nnx + 1 : ii * nnx);
    end
    
    %-----------------------------------------------------------------%
    % NOTE:
    % Up to this point we have computed the matrix containing the
    % approximated solution of the differential problem. However, in
    % order to plot correctly the vector field we must change the point
    % defining the domain as follows. We will create the new matrix X
    % and Y, with the same size, to correctly represent our domain.
    %-----------------------------------------------------------------%
    
    X = zeros(M,nnx);
    Y = zeros(M,nnx);
    horNodes = linspace(cutx(1),cutx(end),nnx);
    inv = 1;

    for kk = 1:nnx

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
    
    %---------------------------------------------------------------------%
    % Create the fitted mesh with quadrature nodes
    %---------------------------------------------------------------------%
    
    M = 250;
    nnx = 250;
    
    Xint = zeros(M,nnx);
    Yint = zeros(M,nnx);
    horNodes = linspace(cutx(1),cutx(end),nnx);
    inv = 1;

    for kk = 1:nnx

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

        Xint(:,kk) = xx;
        Yint(:,kk) = yy;

    end

    
%     % Compute the Gauss-Legendre Nodes
% 
%     objGL = IntegrateHandler();
%     objGL.numbQuadNodes = numbVerInt;
%     [~, refGLNodes, refGLWeights] = gaussLegendre(objGL);
%     
%     % Assign matrix dimensions
%     
%     Xint = zeros(numbVerInt,numbHorInt);
%     Yint = zeros(numbVerInt,numbHorInt);
%     
%     % Create the horizontal nodes by rescalling the reference
%     % Gauss-Legendre nodes
%     
%     a = cutx(1);
%     b = cutx(2);
%     
%     horNodes = a + refGLNodes * (b - a); 
%     
%     % Define parameters
%     
%     inv = 1;
% 
%     for kk = 1:numbHorInt
% 
%         pp  = horNodes(kk);
% 
%         if kk > 1
% 
%             pp1 = horNodes(kk-1);
% 
%             if ((-1/dforma(pp))*(-1/dforma(pp1))<0)
%                 inv = inv + 1;
%             end
% 
%         end
%         
%         % ATTENTION!!!!
%         %-----------------------------------------------------------------%
%         % The original code in the next two lines work only when the point
%         % in the Y direction are symmetric with respect to the X axis. I
%         % perform some modifications to change this, but it still does not
%         % work for the most general case (channel thickness varying).
%         %-----------------------------------------------------------------%  
%         
%         % Create the vertical nodes by rescalling the reference
%         % Gauss-Legendre nodes
%         
%         a1 = pp+(-1)^(inv)*((cuty(2) - cuty(1))/2)*sin(-pi/2+atan(-1/dforma(pp)));
%         b1 = pp-(-1)^(inv)*((cuty(2) - cuty(1))/2)*sin(-pi/2+atan(-1/dforma(pp)));
%      
%         xx = a1 + refGLNodes * (b1 - a1); 
%         xxW = refGLWeights * (b1 - a1); 
%         
%         a2 = ((cuty(2) + cuty(1))/2) + forma(pp)-(-1)^(inv)*((cuty(2) - cuty(1))/2)*cos(-pi/2+atan(-1/dforma(pp)));
%         b2 = ((cuty(2) + cuty(1))/2) + forma(pp)+(-1)^(inv)*((cuty(2) - cuty(1))/2)*cos(-pi/2+atan(-1/dforma(pp)));
%         
%         yy = a2 + refGLNodes * (b2 - a2); 
%         yyW = refGLWeights * (b2 - a2); 
% 
%         Xint(:,kk) = xx;
%         Yint(:,kk) = yy;
%         XintW(:,kk) = xxW;
%         YintW(:,kk) = yyW;
% 
%     end
%     
%     Wint = XintW * YintW;
 
    %---------------------------------------------------------------------%
    % Interpolate and vizualize the freefem++ solution
    %---------------------------------------------------------------------%
    
    [p,e,t,ffSol,stiffMat,massMat] = my_FFimportfilemesh('ffMesh.msh','ffSolution.sol','A.txt','M.txt');
%     ffSolMat = griddata(p(1,:),p(2,:),ffSol,Xint,Yint,'linear');
%  
    higaSol = griddata(X,Y,sol',p(1,:),p(2,:));
    [numbX,numbY] = size(higaSol);
    
    for ii = 1:numbX
        for jj = 1:numbY
            if(isnan(higaSol(ii,jj)))
                higaSol(ii,jj) = 0;
            end
        end
    end
    
%     figure;
%     pdeplot(p,e,t,'XYData',higaSol);
    
    err = higaSol - ffSol;
    
    figure;
    pdeplot(p,e,t,'XYData',higaSol,'ZData',higaSol);
    
    
    figure;
    pdeplot(p,e,t,'XYData',err);
    
    disp('Size Error Vector : ');
    disp(size(err));
    disp('Size Mass Matrix : ');
    disp(size(massMat));
    
    errL2 = sqrt(err * massMat * err');
    errH1 = sqrt(err * massMat * err' + err * stiffMat * err');
    
%     figure;
%     mesh(Xint,Yint,ffSolMat);
%     az = 0;
%     el = 90;
%     view(az, el);
%     title('Freefem++ Solution');
    
    %---------------------------------------------------------------------%
    % Interpolate and vizualize the HigaMod solution
    %---------------------------------------------------------------------%
    
%     interpSol = griddata(X,Y,sol',Xint,Yint);
%     [numbX,numbY] = size(interpSol);
%     
%     for ii = 1:numbX
%         for jj = 1:numbY
%             if(isnan(interpSol(ii,jj)))
%                 interpSol(ii,jj) = 0;
%             end
%         end
%     end
    
    %%%%%%%%%
    % DEBUG %
    %%%%%%%%%
    
%     figure;
%     mesh(Xint,Yint,interpSol);
%     az = 0;
%     el = 90;
%     view(az, el);
%     title('HigaMod Solution');

    %---------------------------------------------------------------------%
    % Import the mass matrix and the stiffness matrix
    %---------------------------------------------------------------------%
    
    
    
    %---------------------------------------------------------------------%
    % Compute error in the mesh points
    %---------------------------------------------------------------------%
    
%     errorMat = interpSol - ffSolMat;
%     
%     [numbX,numbY] = size(errorMat);
%     
%     for ii = 1:numbX
%         for jj = 1:numbY
%             if(isnan(errorMat(ii,jj)))
%                 errorMat(ii,jj) = 0;
%             end
%             if(isnan(ffSolMat(ii,jj)))
%                 ffSolMat(ii,jj) = 0;
%             end
%         end
%     end
%     
%     errL2 = norm(errorMat,2)/norm(ffSolMat,2);
%     errH1 = 0;

    %---------------------------------------------------------------------%
    % Contour plots
    %---------------------------------------------------------------------%
    
%     figure;
%             
%     [~,~] = contourf(Xint,Yint,interpSol,20);    
%     colormap(parula);
%     cmin = min(min(sol'));%-0.05;%min(min(sol'));
%     cmax = max(max(sol')); %0.2;% max(max(sol'));
%     caxis([cmin cmax])
%     colorbar();
%     axis equal
%     set(gca, 'FontSize', 14);
%     
%     figure;
%             
%     [~,~] = contourf(Xint,Yint,ffSolMat,20);    
%     colormap(parula);
%     cmin = min(min(sol'));%-0.05;%min(min(sol'));
%     cmax = max(max(sol')); %0.2;% max(max(sol'));
%     caxis([cmin cmax])
%     colorbar();
%     axis equal
%     set(gca, 'FontSize', 14);
%     
%     figure;
%             
%     [~,~] = contourf(Xint,Yint,errorMat,20);    
%     colormap(parula);
%     cmin = min(min(errorMat));%-0.05;%min(min(sol'));
%     cmax = max(max(errorMat)); %0.2;% max(max(sol'));
%     caxis([cmin cmax])
%     colorbar();
%     axis equal
%     set(gca, 'FontSize', 14);
    
end
