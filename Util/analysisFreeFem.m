function [errorL2,errorH1] = analysisFreeFem(Domain, Mesh,Coeff,force)

    %% Create the Freefem++ simulation folder
    
    for ii = 1:1000
        
        checkFolder = exist(['SimulationResults',num2str(ii)]);
        
        if (checkFolder == 7)
            disp(['The folder SimulationResults',num2str(ii),' already exists!'])
        else
            fileName = ['SimulationResults',num2str(ii)];
            mkdir(fileName)
            break
        end
    end
    
    %% Open the new folder created
    
    cd(fileName);
    
    %% Create file with that will run the Freefem++ simulation
    
    simulationFF = fopen('ffSimulation.edp','w');
    
    %% Write the domain information on the simulation file
    
    fprintf(simulationFF, '// DOMAIN INFORMATION \n\n');
    
    fprintf(simulationFF, ['real minHor =', num2str(Domain.left), ';\n']);
    fprintf(simulationFF, ['real maxHor =', num2str(Domain.right), ';\n']);
    fprintf(simulationFF, ['real minVer =', num2str(Domain.down), ';\n']);
    fprintf(simulationFF, ['real maxVer =', num2str(Domain.up), ';\n\n']);
    
    %% Write the mesh information on the simulation file
    
    fprintf(simulationFF, '// MESH INFORMATION \n\n');
    
    fprintf(simulationFF, ['int numbHorNodes =', num2str(Mesh.numbHorNodes), ';\n']);
    fprintf(simulationFF, ['int numbVerNodes =', num2str(Mesh.numbVerNodes), ';\n\n']);
    
    fprintf(simulationFF, ['real horStep =', num2str((Domain.right - Domain.left)/Mesh.numbHorNodes), ';\n']);
    fprintf(simulationFF, ['real verStep =', num2str((Domain.up - Domain.down)/Mesh.numbVerNodes), ';\n\n']);
    
    %% Write the boundary information on the simulation file
    
    fprintf(simulationFF, '// BOUNDARY INFORMATION \n\n');
    
    fprintf(simulationFF, 'border aa(t=0,1) {x=minHor+t*(maxHor-minHor); y=minVer; label=1;} \n');
    fprintf(simulationFF, 'border bb(t=0,1) {x=maxHor; y=minVer+t*(maxVer-minVer); label=2;} \n');
    fprintf(simulationFF, 'border cc(t=0,1) {x=maxHor-t*(maxHor-minHor); y=maxVer; label=3;} \n');
    fprintf(simulationFF, 'border dd(t=0,1) {x=minHor; y=maxVer-t*(maxVer-minVer); label=4;} \n\n');

    %% Write the mesh creation
    
    fprintf(simulationFF, '// CREATION OF THE MESH \n\n');
    
    fprintf(simulationFF, 'mesh Th=buildmesh( aa(maxHor/horStep) + bb(numbVerNodes) + cc(maxHor/horStep) + dd(numbVerNodes) ); \n\n');
    
    %% Write the bilinear coefficients
    
    fprintf(simulationFF, '// BILINEAR COEFFICIENTS \n\n');
    
    Chi   = func2str(Coeff.chi);
    Mu    = func2str(Coeff.mu);
    Beta1 = func2str(Coeff.beta1);
    Beta2 = func2str(Coeff.beta2);
    Sigma = func2str(Coeff.sigma);
    
    newChi   = strrep(Chi,'@(x,y)',' ');
    newMu    = strrep(Mu,'@(x,y)',' ');
    newBeta1 = strrep(Beta1,'@(x,y)',' ');
    newBeta2 = strrep(Beta2,'@(x,y)',' ');
    newSigma = strrep(Sigma,'@(x,y)',' ');
    
    fprintf(simulationFF, ['func chi =', newChi , ';\n']);
    fprintf(simulationFF, ['func mu =', newMu , ';\n']);
    fprintf(simulationFF, ['func beta1 =', newBeta1 , ';\n']);
    fprintf(simulationFF, ['func beta2 =', newBeta2 , ';\n']);
    fprintf(simulationFF, ['func sigma =', newSigma , ';\n\n']);
    
    %% Write the force term
    
    fprintf(simulationFF, '// DEFINITION OF THE FORCE TERM \n\n');
    
    Force = func2str(force);
    
    newForce = strrep(Force,'@(x,y)',' ');
    newForce = strrep(newForce,'.^','^');
    newForce = strrep(newForce,'.*','*');
    
    fprintf(simulationFF, ['func force =', newForce , ';\n\n']);
    
    %% Write the space definitions and test functions
    
    fprintf(simulationFF, '// SPACE DEFINITION AND TEST FUNCTIONS \n\n');
    
    fprintf(simulationFF, 'fespace Vh(Th,P2); \n');
    fprintf(simulationFF, 'Vh w,u,v,ux,uy; \n');
    fprintf(simulationFF, 'w = force; \n\n');
    
    %% Write the problem definition
    % Note: The boundary condition at the inflow of the channel must be
    % correctly adjusted, because it corrently accepts only homogeneous
    % Dirichlet.
    
    fprintf(simulationFF, '// PROBLEM DEFINITION \n\n');
    
    fprintf(simulationFF, 'problem DirDir(u,v) = \n\n');
    fprintf(simulationFF, 'int2d(Th)     ( mu*( dx(u)*dx(v) + dy(u)*dy(v) ) )   // Laplacian operator \n');
    fprintf(simulationFF, '+ int2d(Th)   ( beta1*dx(u)*v )                      // Advection in x direction \n');
    fprintf(simulationFF, '+ int2d(Th)   ( beta2*dy(u)*v )						// Advection in y direction \n');
    fprintf(simulationFF, '+ int2d(Th)   ( sigma*u*v     )						// Reaction \n');
    fprintf(simulationFF, '- int2d(Th)   ( force*v       )						// Force term \n');
    fprintf(simulationFF, '+ on(4,u=0)									        // Dirichlet Inflow \n');
    fprintf(simulationFF, '+ on(1,u=0)                                          // Dirichlet on Lower Boundary \n');
    fprintf(simulationFF, '+ on(3,u=0);                                         // Dirichlet on Upper Boundary \n\n');
    
    %% Write the command to solve the problem
    
    fprintf(simulationFF, '// SOLUTION OF THE PROBLEM \n\n');
    
    fprintf(simulationFF, 'DirDir; \n\n');
    
    %% Creation of the augmented mesh
    
    fprintf(simulationFF, '// CREATION OF THE AUGMENTED MESH \n\n');
    
    fprintf(simulationFF, 'real[int] meshx(3 * (maxHor - minHor) * numbHorNodes); \n\n');
    fprintf(simulationFF, 'real h = horStep/(maxHor - minHor); \n\n');
    fprintf(simulationFF, 'for (int i=0 ; i < (maxHor - minHor) * numbHorNodes ; i++) \n');
    fprintf(simulationFF, '{ \n');
    fprintf(simulationFF, '    meshx[3*i] = i * h; \n');
    fprintf(simulationFF, '    meshx[3*i+1]=i*h+h/2.; \n');
    fprintf(simulationFF, '    meshx[3*i+2]=i*h+h; \n');
    fprintf(simulationFF, '} \n\n');
    
    %% Write the solution information on the simulation file
    
    fprintf(simulationFF, '// WRITE THE SIMULATION DATA \n\n');
    
    fprintf(simulationFF, '{ \n');
    fprintf(simulationFF, '    ofstream file("DDFF.out"); \n');
    fprintf(simulationFF, '    file.scientific; \n');
    fprintf(simulationFF, '    file.precision(16); \n');
    fprintf(simulationFF, '    int i, j; \n');
    fprintf(simulationFF, '    int M = numbVerNodes; \n');
    fprintf(simulationFF, '    int N = numbHorNodes; \n\n');
    fprintf(simulationFF, '    int NN = 3 * (maxHor - minHor) * N; \n\n');
    fprintf(simulationFF, '    file <<"            "<< NN <<endl; \n');
    fprintf(simulationFF, '    file <<"            "<< numbVerNodes << endl; \n');
    fprintf(simulationFF, '    ux = dx(u); \n');
    fprintf(simulationFF, '    uy = dy(u); \n');
    fprintf(simulationFF, '    for ( int k = 1 ; k <= N * (maxHor - minHor) ; k++) \n');
    fprintf(simulationFF, '    { \n');
    fprintf(simulationFF, '         for (int j = 1 ; j <= numbVerNodes ; j++) \n');
    fprintf(simulationFF, '             file<<"       "<<meshx[ (k-1)*3 + 0]<<"       "<<((j-1.)/(M-1.))*(maxVer-minVer)+minVer<<"       "<<u(meshx[ (k-1)*3 + 0],((j-1.)/(M-1.))*(maxVer-minVer)+minVer)<<"       "<<ux(meshx[ (k-1)*3 + 0],((j-1.)/(M-1.))*(maxVer-minVer)+minVer) <<"       "<<uy(meshx[ (k-1)*3 + 0],((j-1.)/(M-1.))*(maxVer-minVer)+minVer) << endl; \n');
    fprintf(simulationFF, '         for (int j = 1 ; j <= numbVerNodes ; j++) \n');
    fprintf(simulationFF, '         	file<<"       "<<meshx[ (k-1)*3 + 1]<<"       "<<((j-1.)/(M-1.))*(maxVer-minVer)+minVer<<"       "<<u(meshx[ (k-1)*3 + 1],((j-1.)/(M-1.))*(maxVer-minVer)+minVer)<<"       "<<ux(meshx[ (k-1)*3 + 1],((j-1.)/(M-1.))*(maxVer-minVer)+minVer) <<"       "<<uy(meshx[ (k-1)*3 + 1],((j-1.)/(M-1.))*(maxVer-minVer)+minVer) << endl; \n');
    fprintf(simulationFF, '         for (int j = 1 ; j <= numbVerNodes ; j++) \n');
    fprintf(simulationFF, '         	file<<"       "<<meshx[ (k-1)*3 + 2]<<"       "<<((j-1.)/(M-1.))*(maxVer-minVer)+minVer<<"       "<<u(meshx[ (k-1)*3 + 2],((j-1.)/(M-1.))*(maxVer-minVer)+minVer)<<"       "<<ux(meshx[ (k-1)*3 + 2],((j-1.)/(M-1.))*(maxVer-minVer)+minVer) <<"       "<<uy(meshx[ (k-1)*3 + 2],((j-1.)/(M-1.))*(maxVer-minVer)+minVer) << endl; \n');
    fprintf(simulationFF, '    } \n');
    fprintf(simulationFF, '} \n');
    
    %% Close simulation file
    
    fclose(simulationFF);
    
    %%% Create script to find freefem++
    
    %findFF = fopen('findFF','w');
    
    %fprintf(findFF, '#!/bin/bash \n');
    %fprintf(findFF, 'which; \n');
    %fprintf(findFF, 'exit; \n');
    
    %fclose(findFF);
    
    %%% Make the run simulation file executable
    
    %command = 'chmod +x runSimulation';
    %status = system(command);
    
    %if (status == 0)
    %    disp('Simulation is executable.')
    %end
    
    %% Runs the Freefem++ simulation script
    
    command = '/usr/local/ff++/mpich/3.45/bin/freefem++ ffSimulation.edp;';
    status = system(command);
    
    if (status == 0)
        disp('Freefem++ simulation ran without error.')
    else
        disp('Freefem++ simulation stoped.')
    end
    
    %% Copy the simulation analysis file to the current folder
    
    copyfile ../../Util/simulationAnalysis.py
    
    %% Copy the matlab simulation outputs to the current folder
    
    for ii = 1:1000
        
        checkFolder = exist(['MatlabOutput',num2str(ii)]);
        
        if (checkFolder == 7)
            disp(['The folder MatlabOutput',num2str(ii),' is not the last. Skip!'])
        else
            path1 = ['../MatlabOutput',num2str(ii),'/matlab_result.out'];
            path2 = ['../MatlabOutput',num2str(ii),'/quadrature_node.out'];
            
            copyfile(path1);
            copyfile(path2);
            
            break
        end
    end
    
    %% Launch the simulation analysis
    
    command = '/Users/YvesBarbosa/anaconda/bin/python simulationAnalysis.py -inputF . -outputF . -hh matlab_result.out -hq quadrature_node.out -ff DDFF.out -o Result';
    status = system(command);
    status = system(command);
    
    if (status == 0)
        disp('Simulation analysis ran without error.')
    else
        disp('Simulation analysis stoped.')
    end
    
    %% Read analysis report
    
    fileID = fopen('Result.txt','r');
    C = textscan(fileID,'%s');
    fclose(fileID);
    
    strL2 = C{1}(1);
    strH1 = C{1}(2);
    
    [~,errorL2] = strtok(strL2,'=');
    [~,errorH1] = strtok(strH1,'=');
    
    errorL2 = strrep(errorL2,'=','');
    errorH1 = strrep(errorH1,'=','');
    
    errorL2 = str2double(errorL2);
    errorH1 = str2double(errorH1);
    
    %% Exit the simulation folder
    
    cd ..
    
end