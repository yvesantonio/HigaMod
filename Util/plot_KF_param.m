function [errL2,errH1,solsol,meshx,evalNodesY] = plot_KF_param(size_mb,a_ril,b_ril,cutx,cuty,...
                                                                              hx, u,bc_up,bc_down,Coeff_forma,...
                                                                              caso,sol_exx,forma,dforma,psJac,...
                                                                              p,soldx,soldy,k,iteration,timeDomain,...
                                                                              noisySignal,assimilatedState,coeffMean,...
                                                                              coeffVar)

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

    switch caso
            case {1,2,3,4,5,6,7,8,9,10}

            for m = 1:size_mb
                
                approx((m-1)*nnx+1:m*nnx)  = bspeval(p,u((m-1)*nnx+1:m*nnx)',knot,xx3);
                [dc,dk] = bspderiv(p,u((m-1)*nnx+1:m*nnx)',knot);
                dapprox((m-1)*nnx+1:m*nnx) = bspeval(p-1,dc,dk,xx3);
                
                approxN((m-1)*nnx+1:m*nnx)  = bspeval(p,noisySignal((m-1)*nnx+1:m*nnx)',knot,xx3);
                [dc,dk] = bspderiv(p,noisySignal((m-1)*nnx+1:m*nnx)',knot);
                dapproxN((m-1)*nnx+1:m*nnx) = bspeval(p-1,dc,dk,xx3);
                
                approxA((m-1)*nnx+1:m*nnx)  = bspeval(p,assimilatedState((m-1)*nnx+1:m*nnx)',knot,xx3);
                [dc,dk] = bspderiv(p,assimilatedState((m-1)*nnx+1:m*nnx)',knot);
                dapproxA((m-1)*nnx+1:m*nnx) = bspeval(p-1,dc,dk,xx3);

                for j = 2:nnx
                     dapprox((m-1)*nnx+j) = dapprox((m-1)*nnx+j)/sqrt(1+(dforma(xx3(j)))^2);
                     dapproxN((m-1)*nnx+j) = dapproxN((m-1)*nnx+j)/sqrt(1+(dforma(xx3(j)))^2);
                     dapproxA((m-1)*nnx+j) = dapproxA((m-1)*nnx+j)/sqrt(1+(dforma(xx3(j)))^2);
                end

            end;

            % COMPUTATION OF THE NURBS BASE

            u  = approx';
            uN = approxN';
            uA = approxA';

            % COMPUTATION OF THE DERIVATIVE OF THE NURBS BASE

            ux = dapprox';
            uxN = dapproxN';
            uxA = dapproxA';
    end

    solsol = approx';
    solsolN = approxN';
    solsolA = approxA';

    for h = 1:nnx

        for k = 1:M

            for imb = 1:size_mb

                % COMPUTATION OF THE APPROXIMATED SOLUTION VECTOR FIELD
                % EVALUATED IN THE POINTS OF THE DOMAIN MESH

                sol(h,k)   = sol(h,k)   + u(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                
                solN(h,k)   = solN(h,k)   + uN(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                
                solA(h,k)   = solA(h,k)   + uA(h+(imb-1)*nnx)*coeffModalBase(k,imb);

                % COMPUTATION OF THE DERIVATIVE OF THE APPROXIMATED
                % SOLUTION VECTOR FIELD EVALUATED IN THE POINTS OF THE
                % DOMAIN MESH

                sol_x(h,k) = sol_x(h,k) + ux(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                sol_y(h,k) = sol_y(h,k) + u(h+(imb-1)*nnx)*coeffModalBaseDer(k,imb);
                
                sol_xN(h,k) = sol_xN(h,k) + uxN(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                sol_yN(h,k) = sol_yN(h,k) + uN(h+(imb-1)*nnx)*coeffModalBaseDer(k,imb);
                
                sol_xA(h,k) = sol_xA(h,k) + uxA(h+(imb-1)*nnx)*coeffModalBase(k,imb);
                sol_yA(h,k) = sol_yA(h,k) + uA(h+(imb-1)*nnx)*coeffModalBaseDer(k,imb);

            end

        sol(h,k) = sol(h,k)+ a_ril(1) * evalNodesY(k) + b_ril(1);
        
        solN(h,k) = solN(h,k)+ a_ril(1) * evalNodesY(k) + b_ril(1);
        
        solA(h,k) = solA(h,k)+ a_ril(1) * evalNodesY(k) + b_ril(1);

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

%     errL2 = 0;
%     errH1 = 0;
%     aus_figure = 0;
%     solsol = 0;
%     meshx = 0;
%     evalNodesY = 0;
    
    % CREATION OF THE MESH TO EVALUATE THE EXACT SOLUTION OF THE PROBLEM 
        
    sol_ex   = sol_exx(X,Y,timeDomain(iteration));
    sol_exdx = soldx(X,Y,timeDomain(iteration));
    sol_exdy = soldy(X,Y,timeDomain(iteration));
    
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
    
    errL2 = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2)));
    
    % COMPUTATION OF THE H1 NORM OF THE ERROR
    
    errH1 = sqrt(sum(sum(err1.^2))/sum(sum(sol_ex.^2))  +...
                 sum(sum(err2.^2))/sum(sum(sol_exdx.^2))+...
                 sum(sum(err3.^2))/sum(sum(sol_exdy.^2)));
             
    %---------------------------------------------------------------------%
    %       CONTOUR PLOT OF THE APPROXIMATED SOLUTION IN THE DOMAIN
    %---------------------------------------------------------------------%
    
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

    subplot(2,2,1);
    [~,~] = contourf(X,Y,sol',20);    
    colormap(parula);
    cmin = 0;%min(min(sol));
    cmax = 0.15;%max(max(sol));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 10)
    title('Phenomenon')
    
    subplot(2,2,2);
    [~,~] = contourf(X,Y,solN',20);    
    colormap(parula);
    cmin = 0;%min(min(sol));
    cmax = 0.15;%max(max(sol));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 10)
    title('Noisy Measurements')
    
    subplot(2,2,3);
    [~,~] = contourf(X,Y,solA',20);    
    colormap(parula);
    cmin = 0;%min(min(sol));
    cmax = 0.15;%max(max(sol));
    caxis([cmin cmax])
    colorbar();
    axis([minX maxX minY maxY]);
    axis equal
    set(gca, 'FontSize', 10)
    title('Assimilated state')
    
    coeffMeanAug = zeros(size(coeffMean,1) + 1,1);
    coeffMeanAug(2:end) = coeffMean;
    
    grid on
    subplot(2,2,4);
    minX = 1;
    maxX = size(coeffMean,1) + 1;
    minY = 0;
    maxY = 1.2;
    axis([minX maxX minY maxY]);
    axis equal
    plot(ones(length(coeffMeanAug(1:iteration))),'--k')
    plot(coeffMeanAug(1:iteration),'--r*');
    title('Force coefficient')
    
%     subplot(2,1,2);
%     [~,~] = contourf(X,-Y,sol',20);    
%     colormap(parula);
% %     cmin = min(min(sol));
% %     cmax = max(max(sol));
% %     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX -maxY minY]);
%     axis equal
%     set(gca, 'FontSize', 10)
%     title('Lower brush')
    
%     subplot(2,2,3);
%     [~,~] = contourf(X,Y,err1',20);    
%     colormap(parula);
%     cmin = 0;%-0.05;%min(min(sol'));
%     cmax = 0.01;
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 10)
%     title('$\|u_{h} - u\|_{L_{\infty}}$','Interpreter','latex')
%     
%     subplot(2,2,4);
%     [~,~] = contourf(X,Y,err1' + err2' + err3',20);    
%     colormap(parula);
%     cmin = 0;%-0.05;%min(min(sol'));
%     cmax = 0.05;
%     caxis([cmin cmax])
%     colorbar();
%     axis([minX maxX minY maxY]);
%     axis equal
%     set(gca, 'FontSize', 10)
%     title('$\|u_{h} - u\|_{H_{\infty}}$','Interpreter','latex')
%     

    fig2 = figure('visible','off');
    surfSol = surfc(X,Y,solA');    
    az = 55;
    el = 45;
    view(az, el);
    colormap(parula);
%     cmin = 0;%min(min(sol));
%     cmax = 0.15;%max(max(sol));
%     caxis([cmin cmax])
%     colorbar();
    axis off
    axis([0 1 0 1 -0.05 0.25]);
    set(gca, 'FontSize', 10)
    %title('Assimilated state')

    % Save plot

    print(fig1,['Plot_At_it=',num2str(iteration)],'-dpng')
    print(fig2,['Surf_At_it=',num2str(iteration)],'-dpng')
    
    cd ..

end
