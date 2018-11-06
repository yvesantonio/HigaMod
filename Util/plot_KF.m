function [errL2,errH1,meshx,evalNodesY] = plot_KF(size_mb,a_ril,b_ril,cutx,cuty,...
                                                                              hx, u,bc_up,bc_down,Coeff_forma,...
                                                                              caso,sol_exx,forma,dforma,psJac,...
                                                                              p,soldx,soldy,k,iteration,timeDomain,...
                                                                              noisySignal,assimilatedState,Hs,Hm,Hb)

                                                                          
    %% IMPORT CLASSES
            
    import Core.AssemblerADRHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne * k + p + 1 - k;           % Number of Control Point on the X Direction
                                        % for the IGA Mesh

    %% SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 64;
    
    %% DEFINE QUADRATURE PROPERTIES
    
    nx = ne;
    my = M - 1;
    
    qh = 0:nx;
    mh = fliplr(vander(qh));
    for ii = 0: nx
        bh(ii + 1) = nx^(ii + 1)/(ii + 1);
    end
    qh = mh'\bh';
    
    qv = 0:my;
    mv = fliplr(vander(qv));
    for ii = 0: my
        bv(ii + 1) = my^(ii + 1)/(ii + 1);
    end
    qv = mv'\bv';

    [Hgrid,Vgrid] = meshgrid(qh,qv);
    Wmat = Hgrid .* Vgrid;

    %% CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

    evalNodesY = linspace(cuty(1),cuty(2),M);

    y_elements = length(evalNodesY);

    %% INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(nnx,y_elements);
    sol_x = zeros(nnx,y_elements);
    sol_y = zeros(nnx,y_elements);

    solN   = zeros(nnx,y_elements);
    sol_xN = zeros(nnx,y_elements);
    sol_yN = zeros(nnx,y_elements);
    
    solA   = zeros(nnx,y_elements);
    sol_xA = zeros(nnx,y_elements);
    sol_yA = zeros(nnx,y_elements);
    
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

    %% SETTING OF THE ISOGEOMETRIC MESH IN THE X DIRECTION

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

    minX = min(min(min(X),0));
    maxX = max(max(max(X),0));
    minY = min(min(min(Y),0));
    maxY = max(max(max(Y),0));
    
    %% CONTOUR PLOT OF THE EXACT SOLUTION IN THE DOMAIN
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

    %% EVALUATION OF THE SOLUTION ERROR
    
    % EVALUATION OF THE ERROR OF THE SOLUTION AND ITS DERIVATIVE IN THE X
    % AND Y DIRECTIONS
    
    solMat = zeros(M,nnx);
    assMat = zeros(M,nnx);
    
    solR = Hm * Hb * u;
    solA = Hm * Hb * assimilatedState;
    
    for ii = 1:M
        solMat(ii,:) = solR((ii-1)*nnx + 1 : ii*nnx);
        assMat(ii,:) = solA((ii-1)*nnx + 1 : ii*nnx);
    end
    
    numbVerNodes = 128;
    obj_gaussLegendre_2 = IntegrateHandler();
    obj_gaussLegendre_2.numbQuadNodes = numbVerNodes;
    [~, verGLNodes, verWeights] = gaussLegendre(obj_gaussLegendre_2);  
            
    xEval = verGLNodes;
    yEval = verGLNodes;
    
    [XX,YY] = ndgrid(xEval,yEval);
    [XXW,YYW] = ndgrid(verWeights,verWeights);
    WM = XXW .* YYW;
    

    % CREATION OF THE MESH TO EVALUATE THE EXACT SOLUTION OF THE PROBLEM 
        
    ref   = sol_exx(XX,YY,timeDomain(iteration));
    
    higa = griddata(X,Y,solMat,XX,YY,'linear');
    assDat = griddata(X,Y,assMat,XX,YY,'linear');
    
    err1 = abs(ref - higa);
    err2 = abs(ref - assDat);
    
    % COMPUTATION OF THE L2 NORM OF THE ERROR
    
    errL2_higa = sum(sum(WM .* err1));
    errL2_ass =  sum(sum(WM .* err2));
    errL2 = 0;
    errH1 = 0;
        
    disp(['Error L2-Norm Sol.: ',num2str(errL2_higa)])
    disp(['Error L2-Norm Ass.: ',num2str(errL2_ass)])
             
    %% CONTOUR PLOT OF THE APPROXIMATED SOLUTION IN THE DOMAIN
    %---------------------------------------------------------------------%
    
    ref   = sol_exx(X,Y,timeDomain(iteration));
    
    % Create the Freefem++ simulation folder
    
    for ii = 1:1000
        
        checkFolder = exist(['MatlabPlots',num2str(ii)]);
        
        if (checkFolder == 7)
        else
            fileName = ['MatlabPlots',num2str(ii-1)];
            break
        end
    end
    
    % Open the new folder created
    
    cd(fileName);
    
    % Create plot
    
    fig1 = figure('visible','off');
    grid on

    subplot(1,2,1);
    [~,~] = contourf(X,Y,solMat,20);    
    colormap(parula);
    cmin = min(min(solMat));
    cmax = max(max(solMat));
    caxis([cmin cmax]);
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 10)
    %title('Phenomenon')
    
    subplot(1,2,2);
    [~,~] = contourf(X,Y,assMat,20);    
    colormap(parula);
    cmin = min(min(solMat));
    cmax = max(max(solMat));
    caxis([cmin cmax]);
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 10)
    %title('Assimilated state')
        
%     fig1 = figure('visible','off');
%     grid on
% 
%     subplot(1,3,1);
%     [~,~] = contourf(X,Y,solMat,20);    
%     colormap(parula);
%     cmin = min(min(sol));
%     cmax = max(max(sol));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 10)
%     title('Phenomenon')
%     
%     subplot(1,3,2);
%     [~,~] = contourf(X,Y,noiMat,20);    
%     colormap(parula);
%     cmin = min(min(sol));
%     cmax = max(max(sol));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 10)
%     title('Noisy Measurements')
%     
%     subplot(1,3,3);
%     [~,~] = contourf(X,Y,assMat,20);    
%     colormap(parula);
%     cmin = min(min(sol));
%     cmax = max(max(sol));
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 10)
%     title('Assimilated state')

     print(fig1,['Plot_At_it=',num2str(iteration)],'-dpng')
     
%     print(fig2,['Surf_At_it=',num2str(iteration)],'-dpng')
%     
    cd ..

end
