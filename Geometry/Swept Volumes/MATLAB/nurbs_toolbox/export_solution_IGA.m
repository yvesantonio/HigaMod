function [dataFileSurf,errH1,errL2] = export_solution_IGA(size_mb,a_ril,b_ril,cutx,hx,...
                                                                         u,bc_up,bc_down,Coeff_forma,...
                                                                         caso,forma,dforma,p,k,iteration,...
                                                                         psJac,sol_exx,soldx,soldy,psJac2,...
                                                                         timeDomain)

    len = sum(psJac);

    % SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne * k + p + 1 - k;           % Number of Control Point on the X Direction
                                        % for the IGA Mesh

    % SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 40;

    % CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

    evalNodesY = linspace(0,1,M);

    y_elements = length(evalNodesY);

    % INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(2*nnx,y_elements);
    sol_x = zeros(2*nnx,y_elements);
    sol_y = zeros(2*nnx,y_elements);

    import Core.BasisHandler

    obj_newModalBasis = BasisHandler();

    obj_newModalBasis.dimModalBasis = size_mb;
    obj_newModalBasis.evalNodesY = evalNodesY;
    obj_newModalBasis.labelUpBoundCond = bc_up{1};
    obj_newModalBasis.labelDownBoundCond = bc_down{1};
    obj_newModalBasis.coeffForm = Coeff_forma;

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

    in        = 0; 

    switch caso
        case 9
            fi  = 1*len;
        case {1,2,3,4,5,6,7,8,10}
            fi  = 1;
    end

    knot      = augknt([in fi],p+1);
    h         = hx;
    nel       = ne;
    ins       = sort(reshape((h:h:1-h)'*ones(1,k),1,k*(nel-1)));
    ins       = (fi-in)*ins + in;

    [~,knot] = bspkntins(p,in:(fi-in)*1/p:fi,knot,ins);

    %---------------------------------------------------------------------%
    % Note:
    % The vector 'knot' contains the nodes used in the mesh.
    %---------------------------------------------------------------------%

    xx2 = linspace(0,1*len,2*nnx);
    xx3 = linspace(0,1,2*nnx)';

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
        case 9

            for m = 1:size_mb

                approx((m-1)*2*nnx+1:2*m*nnx)  = bspeval(p,u((m-1)*nnx+1:m*nnx)',knot,xx2);
                [dc,dk] = bspderiv(p,u((m-1)*nnx+1:m*nnx)',knot);
                dapprox((m-1)*2*nnx+1:2*m*nnx) = bspeval(p-1,dc,dk,xx2);
                for j = 2:2*nnx
                     dapprox((m-1)*2*nnx+j) = dapprox((m-1)*2*nnx+j);
                end

            end;

            % COMPUTATION OF THE APPROXIMATED SOLUTION

            u  = approx';

            % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
            % SOLUTION

            ux = dapprox';

        case {1,2,3,4,5,6,7,8,10}

            for m = 1:size_mb

                approx((m-1)*2*nnx+1:2*m*nnx)  = bspeval(p,u((m-1)*nnx+1:m*nnx)',knot,xx3);
                [dc,dk] = bspderiv(p,u((m-1)*nnx+1:m*nnx)',knot);
                dapprox((m-1)*2*nnx+1:2*m*nnx) = bspeval(p-1,dc,dk,xx3);

                for j = 2:2*nnx
                     dapprox((m-1)*2*nnx+j) = dapprox((m-1)*2*nnx+j)/sqrt(1+(dforma(xx3(j)))^2);
                end

            end;

            % COMPUTATION OF THE NURBS BASE

            u  = approx';

            % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE

            ux = dapprox';
    end

    for h = 1:2*nnx

        for k = 1:M

            for imb = 1:size_mb

                % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
                % EVALUATED IN THE POINTS OF THE DOMAIN MESH

                sol(h,k)   = sol(h,k)   + u(h+(imb-1)*2*nnx)*coeffModalBase(k,imb);

                % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
                % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
                % DOMAIN MESH

                sol_x(h,k) = sol_x(h,k) + ux(h+(imb-1)*2*nnx)*coeffModalBase(k,imb);
                sol_y(h,k) = sol_y(h,k) + u(h+(imb-1)*2*nnx)*coeffModalBaseDer(k,imb);

            end

        sol(h,k) = sol(h,k)+ a_ril(1) * evalNodesY(k) + b_ril(1);

        end

    end

    %-----------------------------------------------------------------%
    % NOTE:
    % Up to this point we have computed the matrix containing the
    % approximated solution of the differential problem. However, in
    % order to plot correctly the vector field we must change the point
    % defining the domain as follows. We will create the new matrix X
    % and Y, with the same size, to correctly represent our domain.
    %-----------------------------------------------------------------%

    X = zeros(M,2*nnx);
    Y = zeros(M,2*nnx);
    k = linspace(cutx(1),cutx(end),2*nnx);
    inv = 1;

    for kk = 1:2*nnx

        pp  = k(kk);

        if kk > 1

            pp1 = k(kk-1);

            if ((-1/dforma(pp))*(-1/dforma(pp1))<0)
                inv = inv + 1;
            end

        end

        xx = linspace(pp+(-1)^(inv)*0.5*sin(-pi/2+atan(-1/dforma(pp))),pp-(-1)^(inv)*0.5*sin(-pi/2+atan(-1/dforma(pp))),M);
        yy = linspace(forma(pp)-(-1)^(inv)*0.5*cos(-pi/2+atan(-1/dforma(pp))),forma(pp)+(-1)^(inv)*0.5*cos(-pi/2+atan(-1/dforma(pp))),M);

        X(:,kk) = xx;
        Y(:,kk) = yy;

    end
    
    %---------------------------------------------------------------------%
    %                   EXPORT DATA TO VTK PARAVIEW FILE
    %---------------------------------------------------------------------%
    
    % dataFileSurf = 0;
    
        dataFileSurf = sprintf('dataApprox.%d.vtu',iteration - 1);
        TRI = delaunay(X,Y);
        vtktrisurf(TRI,X,Y,sol','ADRConcentration',dataFileSurf);
        
    %---------------------------------------------------------------------%
    %                   EVALUATION OF THE SOLUTION ERROR
    %---------------------------------------------------------------------%
    
    sol_ex      = zeros(M,2*nnx);
    sol_exdx    = zeros(M,2*nnx);
    sol_exdy    = zeros(M,2*nnx);

    % CREATION OF THE MESH TO EVALUATE THE EXACT SOLUTION OF THE PROBLEM 
    
    exactMeshX = zeros(1,2 * nnx);
    exactMeshX(1) = 0;

    for ii = 2 : 2*nnx
        exactMeshX(ii) = exactMeshX(ii - 1) + psJac2(ii - 1);
    end
    
    exactMeshY = linspace(-0.5,0.5,M);

    for ii = 1:M
        
        for j = 1:2*nnx
            
            % SETTING THE POINTS TO EVALUATE THE EXACT SOLUTION
            
            xv          = exactMeshX(j); 
            yv          = exactMeshY(ii);
            
            % EVALUATION OF THE EXACT SOLUTION
            
            sol_ex(ii,j)   = sol_exx(xv,yv,timeDomain(iteration + 1));
            sol_exdx(ii,j) = soldx(xv,yv,timeDomain(iteration + 1));
            sol_exdy(ii,j) = soldy(xv,yv,timeDomain(iteration + 1));
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
