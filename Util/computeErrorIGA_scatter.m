function [errL2,errH1] = computeErrorIGA_scatter(size_mb,a_ril,b_ril,cutx,...
                                                 hx, u,bc_up,bc_down,Coeff_forma,...
                                                 caso,p,k,space,geometry,map)

        % SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne*(p - k) + k + 1;           % Number of Control Point on the X Direction
                                        % for the IGA Mesh

    % SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 2000;

    % CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

    evalNodesY = linspace(0,1,M);

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
    
    xx3 = {linspace(cutx(1),cutx(2),nnx)};
    evalNodesX = linspace(cutx(1),cutx(2),nnx);

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
    
    switch caso
            case {1,2,3,4,5,6,7,8,9,10}
            for m = 1:size_mb
                
                approx((m-1)*nnx+1:m*nnx)  = sp_eval(u((m-1)*nnx+1:m*nnx), space, geometry, xx3);
                dapprox((m-1)*nnx+1:m*nnx) = sp_eval(u((m-1)*nnx+1:m*nnx), space, geometry, xx3,'gradient');

            end

            % COMPUTATION OF THE NURBS BASE

            u  = approx';

            % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE

            ux = dapprox';
    end
    
    solVect = zeros(nnx * M,1);
    solMat = zeros(M,nnx);

    for h = 1:nnx

        for k = 1:M

            for imb = 1:size_mb

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
    
    finalSol = sol;
    for ii = 1:M
        solMat(ii,:) = finalSol((ii-1)*nnx + 1 : ii * nnx);
    end
    
    disp('Finished EXTRACTING THE HIGAMOD SOLUTION');

    X = mapOut(evalNodesX,evalNodesY,map,1);
    Y = mapOut(evalNodesX,evalNodesY,map,2);
 
    %---------------------------------------------------------------------%
    % Interpolate and vizualize the freefem++ solution
    %---------------------------------------------------------------------%
    
    [p,e,t,ffSol,stiffMat,massMat] = my_FFimportfilemesh('ffMesh.msh','ffSolution.sol','A.txt','M.txt');

    disp('Finished EXTRACTING THE FREEFEM++ SOLUTION');
    
    higaSol = griddata(X,Y,solMat,p(1,:),p(2,:));
    [numbX,numbY] = size(higaSol);
    
    for ii = 1:numbX
        for jj = 1:numbY
            if(isnan(higaSol(ii,jj)))
                higaSol(ii,jj) = 0;
            end
        end
    end
    
    err = higaSol - ffSol;
    
    errL2 = sqrt(err * massMat * err');
    errH1 = sqrt(errL2^2 + err * stiffMat * err');
    
    disp('Finished COMPUTING THE ERROR NORMS');
    
    figure;
    pdeplot(p,e,t,'XYData',higaSol,'Contour','on','ColorMap','jet');
    title('HigaMod Solution')
    % export_fig(sprintf('HigaMod'),'-pdf');
    
    figure;
    pdeplot(p,e,t,'XYData',ffSol,'Contour','on','ColorMap','jet');
    title('Freefem++ Solution')
    % export_fig(sprintf('Freefem++'),'-pdf');
    
    figure;
    pdeplot(p,e,t,'XYData',err,'Contour','on','ColorMap','jet');
    title('Solution Error')
    % export_fig(sprintf('Error'),'-pdf');
    
    disp('Finished PLOTTING INFO');
    
end
