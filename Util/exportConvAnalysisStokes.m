function [] = exportConvAnalysisStokes(structError,mVectP,hVectP)

    import Core.AssemblerADRHandler
    import Core.AssemblerStokesHandler
    import Core.AssemblerNSHandler
    import Core.BoundaryConditionHandler
    import Core.IntegrateHandler
    import Core.EvaluationHandler
    import Core.BasisHandler
    import Core.SolverHandler
    
    %% EXTRACT ERROR MATRICES

    matL2ErrorP  = zeros(length(mVectP),length(hVectP));
    matL2ErrorUx = zeros(length(mVectP),length(hVectP));
    matL2ErrorUy = zeros(length(mVectP),length(hVectP));
    matH1ErrorP  = zeros(length(mVectP),length(hVectP));
    matH1ErrorUx = zeros(length(mVectP),length(hVectP));
    matH1ErrorUy = zeros(length(mVectP),length(hVectP));

    for ii = 1:length(mVectP)
        for jj = 1:length(hVectP)

            aux = structError{ii,jj};

            matL2ErrorP(ii,jj)  = aux.errL2_P;
            matL2ErrorUx(ii,jj) = aux.errL2_Ux;
            matL2ErrorUy(ii,jj) = aux.errL2_Uy;
            matH1ErrorP(ii,jj)  = aux.errH1_P;
            matH1ErrorUx(ii,jj) = aux.errH1_Ux;
            matH1ErrorUy(ii,jj) = aux.errH1_Uy;

        end
    end

    %% CREATE CONVERGENCE ANALYSIS PLOTS
    
    %%%%%%%%%%%%%%%%%%%%%%
    % L2 NORM - PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%
    
    figL2P = figure;
    for kk = 1:length(hVectP)
        loglog(mVectP,matL2ErrorP(:,kk),'-o','LineWidth',3,'MarkerSize',3);
        set(gcf, 'Color', 'w');
        set(gca, 'FontSize', 14);
        xlim auto
        ylim auto
        grid on
        hold on
    end

    loglog(mVectP,mVectP.^-1,'--','LineWidth',2);
    loglog(mVectP,mVectP.^-2,'--','LineWidth',2);

    LegendTitles = cell(1,length(hVectP));
    for kk = 1:length(hVectP)
        LegendTitles{kk} = ['h = ' num2str(hVectP(kk))];
    end
    LegendTitles{kk + 1} = ['Order 1'];
    LegendTitles{kk + 2} = ['Order 2'];
    legend(LegendTitles,'Location','northeast')
    set(figL2P, 'Visible', 'on')
    title('Convergence Analysis $\|p - \tilde{p}\|_{L^{2}(\Omega)}$','interpreter','latex')
    cleanfigure('minimumPointsDistance', 0.1);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % L2 NORM - VELOCITY Ux %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    figL2Ux = figure;
    for kk = 1:length(hVectP)
        loglog(mVectP,matL2ErrorUx(:,kk),'-o','LineWidth',3,'MarkerSize',3);
        set(gcf, 'Color', 'w');
        set(gca, 'FontSize', 14);
        xlim auto
        ylim auto
        grid on
        hold on
    end

    loglog(mVectP,mVectP.^-1,'--','LineWidth',2);
    loglog(mVectP,mVectP.^-2,'--','LineWidth',2);

    LegendTitles = cell(1,length(hVectP));
    for kk = 1:length(hVectP)
        LegendTitles{kk} = ['h = ' num2str(hVectP(kk))];
    end
    LegendTitles{kk + 1} = ['Order 1'];
    LegendTitles{kk + 2} = ['Order 2'];
    legend(LegendTitles,'Location','northeast')
    set(figL2Ux, 'Visible', 'on')
    title('Convergence Analysis $\|u_{x} - \tilde{u}_{x}\|_{L^{2}(\Omega)}$','interpreter','latex')
    cleanfigure('minimumPointsDistance', 0.1);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % L2 NORM - VELOCITY Uy %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    figL2Uy = figure;
    for kk = 1:length(hVectP)
        loglog(mVectP,matL2ErrorUy(:,kk),'-o','LineWidth',3,'MarkerSize',3);
        set(gcf, 'Color', 'w');
        set(gca, 'FontSize', 14);
        xlim auto
        ylim auto
        grid on
        hold on
    end

    loglog(mVectP,mVectP.^-1,'--','LineWidth',2);
    loglog(mVectP,mVectP.^-2,'--','LineWidth',2);

    LegendTitles = cell(1,length(hVectP));
    for kk = 1:length(hVectP)
        LegendTitles{kk} = ['h = ' num2str(hVectP(kk))];
    end
    LegendTitles{kk + 1} = ['Order 1'];
    LegendTitles{kk + 2} = ['Order 2'];
    legend(LegendTitles,'Location','northeast')
    set(figL2Uy, 'Visible', 'on')
    title('Convergence Analysis $\|u_{y} - \tilde{u}_{y}\|_{L^{2}(\Omega)}$','interpreter','latex')
    cleanfigure('minimumPointsDistance', 0.1);

    %%%%%%%%%%%%%%%%%%%%%%
    % H1 NORM - PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%

    figH1P = figure;
    for kk = 1:length(hVectP)
        loglog(mVectP,matH1ErrorP(:,kk),'-o','LineWidth',3,'MarkerSize',3);
        set(gcf, 'Color', 'w');
        set(gca, 'FontSize', 14);
        xlim auto
        ylim auto
        grid on
        hold on
    end

    loglog(mVectP,mVectP.^-1,'--','LineWidth',2);
    loglog(mVectP,mVectP.^-2,'--','LineWidth',2);

    LegendTitles = cell(1,length(hVectP));
    for kk = 1:length(hVectP)
        LegendTitles{kk} = ['h = ' num2str(hVectP(kk))];
    end
    LegendTitles{kk + 1} = ['Order 1'];
    LegendTitles{kk + 2} = ['Order 2'];
    legend(LegendTitles,'Location','northeast')
    set(figH1P, 'Visible', 'on')
    title('Convergence Analysis $\|p - \tilde{p}\|_{H^{1}(\Omega)}$','interpreter','latex')
    cleanfigure('minimumPointsDistance', 0.1);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % H1 NORM - VELOCITY Ux %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    figH1Ux = figure;
    for kk = 1:length(hVectP)
        loglog(mVectP,matL2ErrorUx(:,kk),'-o','LineWidth',3,'MarkerSize',3);
        set(gcf, 'Color', 'w');
        set(gca, 'FontSize', 14);
        xlim auto
        ylim auto
        grid on
        hold on
    end

    loglog(mVectP,mVectP.^-1,'--','LineWidth',2);
    loglog(mVectP,mVectP.^-2,'--','LineWidth',2);

    LegendTitles = cell(1,length(hVectP));
    for kk = 1:length(hVectP)
        LegendTitles{kk} = ['h = ' num2str(hVectP(kk))];
    end
    LegendTitles{kk + 1} = ['Order 1'];
    LegendTitles{kk + 2} = ['Order 2'];
    legend(LegendTitles,'Location','northeast')
    set(figH1Ux, 'Visible', 'on')
    title('Convergence Analysis $\|u_{x} - \tilde{u}_{x}\|_{H^{1}(\Omega)}$','interpreter','latex')
    cleanfigure('minimumPointsDistance', 0.1);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % H1 NORM - VELOCITY Uy %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    figH1Uy = figure;
    for kk = 1:length(hVectP)
        loglog(mVectP,matH1ErrorUy(:,kk),'-o','LineWidth',3,'MarkerSize',3);
        set(gcf, 'Color', 'w');
        set(gca, 'FontSize', 14);
        xlim auto
        ylim auto
        grid on
        hold on
    end

    loglog(mVectP,mVectP.^-1,'--','LineWidth',2);
    loglog(mVectP,mVectP.^-2,'--','LineWidth',2);

    LegendTitles = cell(1,length(hVectP));
    for kk = 1:length(hVectP)
        LegendTitles{kk} = ['h = ' num2str(hVectP(kk))];
    end
    LegendTitles{kk + 1} = ['Order 1'];
    LegendTitles{kk + 2} = ['Order 2'];
    legend(LegendTitles,'Location','northeast')
    set(figH1Uy, 'Visible', 'on')
    title('Convergence Analysis $\|u_{y} - \tilde{u}_{y}\|_{H^{1}(\Omega)}$','interpreter','latex')
    cleanfigure('minimumPointsDistance', 0.1);

    %% EXPORT PLOTS

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CREATE NEW EXPORT FOLDER %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd('Results')

    for ii = 1:100        
        checkFolder = exist(['ConvergenceAnalysis',num2str(ii)]);

        if (checkFolder ~= 7)
            fileNameF = ['ConvergenceAnalysis',num2str(ii)];
            mkdir(fileNameF)
            break
        end
    end

    cd(fileNameF)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE L2 NORM PRESSURE RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print(figL2P,'ConvAnalysis_L2_P','-dpng');
    matlab2tikz('ConvAnalysis_L2_P.tex','width', '\fwidth','standalone', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE L2 NORM VELOCITY Ux RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print(figL2Ux,'ConvAnalysis_L2_Ux','-dpng');
    matlab2tikz('ConvAnalysis_L2_Ux.tex','width', '\fwidth','standalone', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE L2 NORM VELOCITY Uy RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print(figL2Uy,'ConvAnalysis_L2_Uy','-dpng');
    matlab2tikz('ConvAnalysis_L2_Uy.tex','width', '\fwidth','standalone', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE H1 NORM PRESSURE RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print(figH1P,'ConvAnalysis_H1_P','-dpng');
    matlab2tikz('ConvAnalysis_H1_P.tex','width', '\fwidth','standalone', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE H1 NORM VELOCITY Ux RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print(figH1Ux,'ConvAnalysis_H1_Ux','-dpng');
    matlab2tikz('ConvAnalysis_H1_Ux.tex','width', '\fwidth','standalone', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE H1 NORM VELOCITY Uy RESULTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print(figH1Uy,'ConvAnalysis_H1_Uy','-dpng');
    matlab2tikz('ConvAnalysis_H1_Uy.tex','width', '\fwidth','standalone', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % CLOSE CURRENT FOLDER %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    cd ../..

end
