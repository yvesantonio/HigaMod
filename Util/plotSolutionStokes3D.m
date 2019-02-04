function [] = plotSolutionStokes(plotStruct)

    import Core.AssemblerADRHandler
    import Core.AssemblerStokesHandler
    import Core.AssemblerNSHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINTION HORIZONTAL NODES AND INTERVALS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numbKnotsP      = plotStruct.numbKnotsP;
    numbKnotsUx     = plotStruct.numbKnotsUx;
    numbKnotsUy     = plotStruct.numbKnotsUy;
    numbKnotsUz     = plotStruct.numbKnotsUz;
    
    numbControlPtsP  = plotStruct.numbControlPtsP;
    numbControlPtsUx = plotStruct.numbControlPtsUx;
    numbControlPtsUy = plotStruct.numbControlPtsUy;
    numbControlPtsUz = plotStruct.numbControlPtsUz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINITION NUMBER OF VERTICAL EVALUATION POINTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M = 50;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINITION VERTICAL EVALUATION POINTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    evalNodesY     = linspace(0,1,M);
    numbEvalNodesY = length(evalNodesY);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIALIZATION OF THE SOLUTION MATRICES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    solP   = zeros(numbControlPtsP ,numbEvalNodesY,numbEvalNodesY);
    solUx  = zeros(numbControlPtsUx,numbEvalNodesY,numbEvalNodesY);
    solUy  = zeros(numbControlPtsUy,numbEvalNodesY,numbEvalNodesY);
    solUz  = zeros(numbControlPtsUz,numbEvalNodesY,numbEvalNodesY);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % PRESSURE MODAL BASIS %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objNewModalBasisP = BasisHandler();

    objNewModalBasisP.dimModalBasis         = plotStruct.discStruct.numbModesP;
    objNewModalBasisP.evalNodesY            = evalNodesY;
    objNewModalBasisP.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_P;
    objNewModalBasisP.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_P;
    objNewModalBasisP.coeffForm             = plotStruct.probParameters;

    [modalBasisP,~,~] = newModalBasisStokes3D(objNewModalBasisP);
    
    %%%%%%%%%%%%%%%%%%
    % Ux MODAL BASIS %
    %%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objNewModalBasisUx = BasisHandler();

    objNewModalBasisUx.dimModalBasis         = plotStruct.discStruct.numbModesUx;
    objNewModalBasisUx.evalNodesY            = evalNodesY;
    objNewModalBasisUx.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_Ux;
    objNewModalBasisUx.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_Ux;
    objNewModalBasisUx.coeffForm             = plotStruct.probParameters;

    [modalBasisUx,~,~] = newModalBasisStokes3D(objNewModalBasisUx);
    
    %%%%%%%%%%%%%%%%%%
    % Uy MODAL BASIS %
    %%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objNewModalBasisUy = BasisHandler();

    objNewModalBasisUy.dimModalBasis         = plotStruct.discStruct.numbModesUy;
    objNewModalBasisUy.evalNodesY            = evalNodesY;
    objNewModalBasisUy.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_Uy;
    objNewModalBasisUy.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_Uy;
    objNewModalBasisUy.coeffForm             = plotStruct.probParameters;

    [modalBasisUy,~,~] = newModalBasisStokes3D(objNewModalBasisUy);
    
    %%%%%%%%%%%%%%%%%%
    % Uz MODAL BASIS %
    %%%%%%%%%%%%%%%%%%
    
    import Core.BasisHandler

    objNewModalBasisUz = BasisHandler();

    objNewModalBasisUz.dimModalBasis         = plotStruct.discStruct.numbModesUz;
    objNewModalBasisUz.evalNodesY            = evalNodesY;
    objNewModalBasisUz.labelUpBoundCond      = plotStruct.boundaryStruct.boundCondStruc.bc_up_tag_Uz;
    objNewModalBasisUz.labelDownBoundCond    = plotStruct.boundaryStruct.boundCondStruc.bc_down_tag_Uz;
    objNewModalBasisUz.coeffForm             = plotStruct.probParameters;

    [modalBasisUz,~,~] = newModalBasisStokes3D(objNewModalBasisUz);


    %---------------------------------------------------------------------%
    % Note:
    % The coefficients 'coeffModalBase' and 'coeffModalBaseDer' are
    % vectors conrresponding respectively to the coefficients of the
    % modal bases and derivative of the modal basis on the points
    % assigned in the Y direction.
    %---------------------------------------------------------------------%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRESSURE HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXP = {linspace(0,1,numbControlPtsP)};
    evalNodesXP       = linspace(0,1,numbControlPtsP);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Ux HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXUx = {linspace(0,1,numbControlPtsUx)};
    evalNodesXUx       = linspace(0,1,numbControlPtsUx);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Uy HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXUy = {linspace(0,1,numbControlPtsUy)};
    evalNodesXUy       = linspace(0,1,numbControlPtsUy);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Uz HORIZONTAL MESH %
    %%%%%%%%%%%%%%%%%%%%%%
    
    structEvalNodesXUz = {linspace(0,1,numbControlPtsUz)};
    evalNodesXUz       = linspace(0,1,numbControlPtsUz);

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTRACT SOLUTION VECTOR %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    AUX = plotStruct.solStruct;

    index1 = numbControlPtsUx * plotStruct.discStruct.numbModesUx;
    index2 = numbControlPtsUy * plotStruct.discStruct.numbModesUy;
    index3 = numbControlPtsUz * plotStruct.discStruct.numbModesUz;
    index4 = numbControlPtsP  * plotStruct.discStruct.numbModesP;

    uUx = AUX(1:index1,1);
    uUy = AUX(index1 + 1 : index1 + index2,1);
    uUz = AUX(index1 + index2 + 1 : index1 + index2 + index3,1);
    uP  = AUX(index1 + index2 + index3 + 1 : index1 + index2 + index3 + index4,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % SOLUTION PROJECTION %
    %%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Project the solution vector of the pressure using the
    % information of the isogeometric basis and the modal basis.
    %---------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COORDINATES OF THE PROJECTION IN THE PHYSICAL DOMAIN %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PROJECTION POINTS
    
    [XP ,YP ,ZP ] = mapOut3DHiMod(evalNodesXP ,evalNodesY,evalNodesY,plotStruct.geometricInfo,plotStruct.geometricInfo.Type);
    [XUx,YUx,ZUx] = mapOut3DHiMod(evalNodesXUx,evalNodesY,evalNodesY,plotStruct.geometricInfo,plotStruct.geometricInfo.Type);
    [XUy,YUy,ZUy] = mapOut3DHiMod(evalNodesXUy,evalNodesY,evalNodesY,plotStruct.geometricInfo,plotStruct.geometricInfo.Type);
    [XUz,YUz,ZUz] = mapOut3DHiMod(evalNodesXUz,evalNodesY,evalNodesY,plotStruct.geometricInfo,plotStruct.geometricInfo.Type);

    % FINE MESH FOR THE VISUALIZATION OF THE SOLUTION

    Nx = 10;
    Ny = 10;
    Nz = 10;
    
    x_eval = linspace(0,1,Nx);
    y_eval = linspace(0,1,Ny);
    z_eval = linspace(0,1,Nz);
    
    [Xeval,Yeval,Zeval] = mapOut3DHiMod(x_eval,y_eval,z_eval,plotStruct.geometricInfo,plotStruct.geometricInfo.Type);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROJECTION USING THE ISOGEOMETRIC BASIS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for m = 1:plotStruct.discStruct.numbModesP

        currUP = uP;
        spaceP = plotStruct.spaceP;
        geoP   = plotStruct.refDomain1D;
        approxP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP)  = sp_eval(currUP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP), spaceP, geoP, structEvalNodesXP);
        dapproxP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP) = sp_eval(currUP((m-1) * numbControlPtsP + 1 : m * numbControlPtsP), spaceP, geoP, structEvalNodesXP,'gradient');

    end

    for m = 1:plotStruct.discStruct.numbModesUx

        currUUx = uUx;
        spaceUx = plotStruct.spaceUx;
        geoUx   = plotStruct.refDomain1D;
        approxUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx)  = sp_eval(currUUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx), spaceUx, geoUx, structEvalNodesXUx);
        dapproxUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx) = sp_eval(currUUx((m-1) * numbControlPtsUx + 1 : m * numbControlPtsUx), spaceUx, geoUx, structEvalNodesXUx,'gradient');

    end

    for m = 1:plotStruct.discStruct.numbModesUy

        currUUy = uUy;
        spaceUy = plotStruct.spaceUy;
        geoUy   = plotStruct.refDomain1D;
        approxUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy)  = sp_eval(currUUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy), spaceUy, geoUy, structEvalNodesXUy);
        dapproxUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy) = sp_eval(currUUy((m-1) * numbControlPtsUy + 1 : m * numbControlPtsUy), spaceUy, geoUy, structEvalNodesXUy,'gradient');

    end
    
    for m = 1:plotStruct.discStruct.numbModesUz

        currUUz = uUz;
        spaceUz = plotStruct.spaceUz;
        geoUz   = plotStruct.refDomain1D;
        approxUz((m-1) * numbControlPtsUz + 1 : m * numbControlPtsUz)  = sp_eval(currUUz((m-1) * numbControlPtsUz + 1 : m * numbControlPtsUz), spaceUz, geoUz, structEvalNodesXUz);
        dapproxUz((m-1) * numbControlPtsUz + 1 : m * numbControlPtsUz) = sp_eval(currUUz((m-1) * numbControlPtsUz + 1 : m * numbControlPtsUz), spaceUz, geoUz, structEvalNodesXUz,'gradient');

    end

    igaProjP  = approxP';
    igaProjUx = approxUx';
    igaProjUy = approxUy';
    igaProjUz = approxUz';

    igaProjP_dx  = dapproxP';
    igaProjUx_dx = dapproxUx';
    igaProjUy_dx = dapproxUy';
    igaProjUz_dx = dapproxUz';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATRIX STRUCTURE FOR THE PROJECTED INFORMATION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solMatP   = zeros(M,M,numbControlPtsP);
    solMatUx  = zeros(M,M,numbControlPtsUx);
    solMatUy  = zeros(M,M,numbControlPtsUy);
    solMatUz  = zeros(M,M,numbControlPtsUz);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROJECTION USING THE MODAL BASIS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for kk = 1:M
        for jj = 1:M

            % PROJECTION OF THE PRESSURE

            for hh = 1:numbControlPtsP
                for mb = 1:plotStruct.discStruct.numbModesP
                    solP(hh,kk,jj) = solP(hh,kk,jj) + igaProjP(hh + (mb-1) * numbControlPtsP) * modalBasisP(kk,jj,mb);
                end
            end

            % PROJECTION OF Ux

            for hh = 1:numbControlPtsUx
                for mb = 1:plotStruct.discStruct.numbModesUx
                    solUx(hh,kk,jj) = solUx(hh,kk,jj) + igaProjUx(hh + (mb-1) * numbControlPtsUx) * modalBasisUx(kk,jj,mb);
                end
            end

            % PROJECTION OF Uy

            for hh = 1:numbControlPtsUy
                for mb = 1:plotStruct.discStruct.numbModesUy
                    solUy(hh,kk,jj) = solUy(hh,kk,jj) + igaProjUy(hh + (mb-1) * numbControlPtsUy) * modalBasisUy(kk,jj,mb);
                end
            end
            
            % PROJECTION OF Uz

            for hh = 1:numbControlPtsUz
                for mb = 1:plotStruct.discStruct.numbModesUz
                    solUz(hh,kk,jj) = solUz(hh,kk,jj) + igaProjUz(hh + (mb-1) * numbControlPtsUz) * modalBasisUz(kk,jj,mb);
                end
            end
            
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE THE SOLUTION MATRIX %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    finalSolP  = solP;
    finalSolUx = solUx;
    finalSolUy = solUy;
    finalSolUz = solUz;
    
    for ii = 1:M
        for jj = 1:M
            
            % PRESSURE
            
            for kk = 1:numbControlPtsP
                solMatP(ii,jj,kk) = finalSolP(kk,ii,jj);
            end
            
            % VELOCITY IN X
            
            for kk = 1:numbControlPtsUx
                solMatUx(ii,jj,kk) = finalSolUx(kk,ii,jj);
            end
            
            % VELOCITY IN Y
            
            for kk = 1:numbControlPtsUy
                solMatUy(ii,jj,kk) = finalSolUy(kk,ii,jj);
            end
            
            % VELOCITY IN Z
            
            for kk = 1:numbControlPtsUz
                solMatUz(ii,jj,kk) = finalSolUz(kk,ii,jj);
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROJECTION IN FINE MESH FOR THE PLOTTING %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : This step is necessary because during the previous projection
    % step we could not choose the number of points to plot the solution
    % along the centerline, only along the transverse direction.
    %---------------------------------------------------------------------%
    
    higaSolP  = griddata(XP ,YP ,ZP ,solMatP ,Xeval,Yeval,Zeval);
    higaSolUx = griddata(XUx,YUx,ZUx,solMatUx,Xeval,Yeval,Zeval);
    higaSolUy = griddata(XUy,YUy,ZUy,solMatUy,Xeval,Yeval,Zeval);
    higaSolUz = griddata(XUz,YUz,ZUz,solMatUz,Xeval,Yeval,Zeval);
    
    [numbX,numbY,numbZ] = size(higaSolP);
    
    for ii = 1:numbX
        for jj = 1:numbY
            for kk = 1:numbZ
                
                % PRESSURE
                
                if(isnan(higaSolP(ii,jj,kk)))
                    higaSolP(ii,jj,kk) = 0;
                end
                
                % VELOCITY IN Ux
                
                if(isnan(higaSolUx(ii,jj,kk)))
                    higaSolUx(ii,jj,kk) = 0;
                end
                
                % VELOCITY IN Uy
                
                if(isnan(higaSolUy(ii,jj,kk)))
                    higaSolUy(ii,jj,kk) = 0;
                end
                
                % VELOCITY IN Uz
                
                if(isnan(higaSolUz(ii,jj,kk)))
                    higaSolUz(ii,jj,kk) = 0;
                end
                
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STORE SOLUTION STRUCTURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solStruct.P  = higaSolP;
    solStruct.Ux = higaSolUx;
    solStruct.Uy = higaSolUy;
    solStruct.Uz = higaSolUz;
    solStruct.X  = Xeval;
    solStruct.Y  = Yeval;
    solStruct.Z  = Zeval;

    %%%%%%%%%%%%%%%%%%%%%%%
    % CHECK EXPORT FOLDER %
    %%%%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Create the export folder to save the plots and video of the
    % simulation.
    %---------------------------------------------------------------------%
    
    for ii = 1:1000

        checkFolder = exist(['MatlabPlots',num2str(ii)]);

        if (checkFolder == 7)
        else
            fileName = ['MatlabPlots',num2str(ii-1)];
            break
        end
    end

    %%%%%%%%%%%%%%%%%%%%%
    % OPEN CURRENT FOLDER 
    %%%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Open the current export folder.
    %---------------------------------------------------------------------%
    
    cd(fileName);

    %%%%%%%%%%%%%%%%%%%%
    % CREATE FIGURE PLOT
    %%%%%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Create figure to store the plots and exporting format.
    %---------------------------------------------------------------------%
        
    X  = solStruct.X;
    Y  = solStruct.Y;
    Z  = solStruct.Z;
    P  = solStruct.P;
    Ux = solStruct.Ux;
    Uy = solStruct.Uy;
    Uz = solStruct.Uz;
    U  = sqrt(Ux.^2 + Uy.^2 + Uz.^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SCATTER PLOT OF THE MAGNITUDE OF THE VELOCITY FIELD %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fig1 = figure('visible','on');
    scatter3(XP(:),YP(:),ZP(:),40,solMatP(:),'filled');
    
    colormap(jet);
    
%     cmin = min(min(solMatP(:)));
%     cmax = max(max(solMatP(:)));
%     caxis([cmin cmax]);
    cb = colorbar();
    cb.Label.String = 'Velocity Magnitude';
    
    minX = min(min(X));
    maxX = max(max(X));
    minY = min(min(Y));
    maxY = max(max(Y));
    minZ = min(min(Z));
    maxZ = max(max(Z));
    
    view(45,45)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    set(gca, 'FontSize', 14);
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SCATTER PLOT OF THE MAGNITUDE OF THE VELOCITY FIELD %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fig2 = figure('visible','on');
    scatter3(X(:),Y(:),Z(:),40,Uy(:),'filled');
    
    colormap(jet);
    
    cmin = min(min(Uy(:)));
    cmax = max(max(Uy(:)));
    caxis([cmin cmax]);
    cb = colorbar();
    cb.Label.String = 'Velocity Magnitude';
    
    minX = min(min(X));
    maxX = max(max(X));
    minY = min(min(Y));
    maxY = max(max(Y));
    minZ = min(min(Z));
    maxZ = max(max(Z));
    
    view(45,45)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    set(gca, 'FontSize', 14);
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SCATTER PLOT OF THE MAGNITUDE OF THE VELOCITY FIELD %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fig3 = figure('visible','on');
    scatter3(X(:),Y(:),Z(:),40,Uz(:),'filled');
    
    colormap(jet);
    
    cmin = min(min(Uz(:)));
    cmax = max(max(Uz(:)));
    caxis([cmin cmax]);
    cb = colorbar();
    cb.Label.String = 'Velocity Magnitude';
    
    minX = min(min(X));
    maxX = max(max(X));
    minY = min(min(Y));
    maxY = max(max(Y));
    minZ = min(min(Z));
    maxZ = max(max(Z));
    
    view(45,45)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    set(gca, 'FontSize', 14);
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SCATTER PLOT OF THE MAGNITUDE OF THE PRESSURE FIELD %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     fig2 = figure('visible','off');
%     scatter3(X,Y,Z,40,P,'filled');
%     
%     colormap(jet);
%     
%     cmin = min(min(P));
%     cmax = max(max(P));
%     caxis([cmin cmax]);
%     cb = colorbar();
%     cb.Label.String = 'Velocity Magnitude';
%     
%     minX = min(min(X));
%     maxX = max(max(X));
%     minY = min(min(Y));
%     maxY = max(max(Y));
%     minZ = min(min(Z));
%     maxZ = max(max(Z));
%     
%     view(45,45)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     axis([minX maxX minY maxY minZ maxZ]);
%     axis equal
%     scale = 0.1;
%     daspect([1 1 scale]);
%     pbaspect([1 1 scale]);
%     set(gca, 'FontSize', 14);
%     hold on;

    % EXPORT VTK FILE OF THE HIGAMOD SOLUTION
    
    filename = 'HigaModSol.vtk';
    title    = 'Velocity';
    title2   = 'Pressure';
    
    vtkwrite(filename,'unstructured_grid',X,Y,Z,'vectors',title,Ux,Uy,Uz,'scalars',title2,P);

    %%%%%%%%%%%%%
    % SAVE FIGURE
    %%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : Save the current figure.
    %---------------------------------------------------------------------%

%     print(fig1,'VelocityField','-dpng');
%     print(fig2,'PressureField','-dpng');
    
    cd ..
    
    %%%%%%%%%%%%%%%%%
    % COMPUTE ERROR %
    %%%%%%%%%%%%%%%%%
    %---------------------------------------------------------------------%
    % Note : This function only plots the HigaMod solution and its
    % gradient. We do not compute the error with respect to a analytical or
    % reference solution, so the L2 and H1 norms of the error are set to
    % zero.
    %---------------------------------------------------------------------%

end
