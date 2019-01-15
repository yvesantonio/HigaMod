function [structffP,structffE,structffT,structffSolP,structffSolUx,structffSolUy,structffAp,structffMp,structffMesh] = importSolutionStokesUnsteady(file_Mesh,file_ffSolP,file_ffSolUx,file_ffSolUy,file_Ap,file_Mp,timeDomain)

structffP     = cell(1,length(timeDomain));
structffE     = cell(1,length(timeDomain));
structffT     = cell(1,length(timeDomain));
structffSolP  = cell(1,length(timeDomain));
structffSolUx = cell(1,length(timeDomain));
structffSolUy = cell(1,length(timeDomain));
structffAp    = cell(1,length(timeDomain));
structffMp    = cell(1,length(timeDomain));
structffMesh  = cell(1,length(timeDomain));

for tt = 1:length(timeDomain)
    
    disp(' ')
    disp('-----------------------------------------------------------------------------------------')
    disp(['  START EXTRACTING THE FREEFEM++ SOLTUION FOR THE STOKES PROBLEM FOR K = ',num2str(tt),' '])
    disp('-----------------------------------------------------------------------------------------')
    disp(' ')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IMPORT MESH INFORMATION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (tt == 1 )
        tic;

        rawData1 = importdata(file_Mesh{tt});

        np=rawData1(1,1);
        p = rawData1(2:np+1,1:2)';

        nt=rawData1(1,2);
        t= [rawData1(np+2:2:np+1+2*nt,:)'; rawData1(np+3:2:np+2+2*nt,1)'];

        ne = rawData1(1,3);
        e = [rawData1(np+2+2*nt:np+1+2*nt+ne,1:2)'; zeros(1,ne); ones(1,ne); ...
            rawData1(np+2+2*nt:np+1+2*nt+ne,3)'; ones(1,ne); zeros(1,ne)];

        structMesh = [];
        structMesh.NumbNodes = np;
        structMesh.NumbElements = ne;
        structMesh.Points = p;
        structMesh.Elements = e;
        structMesh.ElementOrder = 4;

        structffP{tt}    = p;
        structffE{tt}    = e;
        structffT{tt}    = t;
        structffMesh{tt} = structMesh;

        t1 = toc;

        disp(['MESH    - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t1),' [s]']);
    end

    %%%%%%%%%%%%%%%%%%%%%
    % IMPORT SOLUTION P %
    %%%%%%%%%%%%%%%%%%%%%

    tic;

    rawData2 = importdata(file_ffSolP{tt});

    solStrVect = rawData2.textdata(6:(end-1));
    temp = strjoin(solStrVect);
    ffSolP = str2num(temp);
    structffSolP{tt} = ffSolP;

    t2 = toc;

    disp(['SOL. P  - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t2),' [s]']);

    %%%%%%%%%%%%%%%%%%%%%%
    % IMPORT SOLUTION Ux %
    %%%%%%%%%%%%%%%%%%%%%%

    tic;

    rawData3 = importdata(file_ffSolUx{tt});

    solStrVect = rawData3.textdata(6:(end-1));
    temp = strjoin(solStrVect);
    ffSolUx = str2num(temp);
    structffSolUx{tt} = ffSolUx;

    t3 = toc;

    disp(['SOL. Ux - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t3),' [s]']);

    %%%%%%%%%%%%%%%%%%%%%%
    % IMPORT SOLUTION Uy %
    %%%%%%%%%%%%%%%%%%%%%%

    tic;

    rawData4 = importdata(file_ffSolUy{tt});

    solStrVect = rawData4.textdata(6:(end-1));
    temp = strjoin(solStrVect);
    ffSolUy = str2num(temp);
    structffSolUy{tt} = ffSolUy;

    t4 = toc;

    disp(['SOL. Uy - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t4),' [s]']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IMPORT STIFFNESS MATRIX FOR PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (tt == 1)
        tic;

        rawData5 = importdata(file_Ap{tt});

        rows    = rawData5.data(2:end,1);
        columns = rawData5.data(2:end,2);
        values  = rawData5.data(2:end,3);

        Ap = sparse(rows,columns,values);
        structffAp{tt} = Ap;

        t5 = toc;

        disp(['STF. P  - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t5),' [s]']);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IMPORT MASS MATRIX FOR PRESSURE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (tt == 1)
        tic;

        rawData6 = importdata(file_Mp{tt});

        rows    = rawData6.data(2:end,1);
        columns = rawData6.data(2:end,2);
        values  = rawData6.data(2:end,3);

        Mp = sparse(rows,columns,values);
        structffMp{tt} = Mp;

        t6 = toc;

        disp(['MSS. P  - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t6),' [s]']);
    end
    
end

end


