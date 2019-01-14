function [p,e,t,ffSolP,ffSolUx,ffSolUy,Ap,Mp,structMesh] = importSolutionStokes(file_Mesh,file_ffSolP,file_ffSolUx,file_ffSolUy,file_Ap,file_Mp,file_Au,file_Mu)

disp(' ')
disp('------------------------------------------------------------------')
disp('  START EXTRACTING THE FREEFEM++ SOLTUION FOR THE STOKES PROBLEM  ')
disp('------------------------------------------------------------------')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT MESH INFORMATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

rawData1 = importdata(file_Mesh);

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

t1 = toc;

disp(['MESH    - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t1),' [s]']);

%%%%%%%%%%%%%%%%%%%%%
% IMPORT SOLUTION P %
%%%%%%%%%%%%%%%%%%%%%

tic;

rawData2 = importdata(file_ffSolP);

solStrVect = rawData2.textdata(6:(end-1));
temp = strjoin(solStrVect);
ffSolP = str2num(temp);

t2 = toc;

disp(['SOL. P  - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t2),' [s]']);

%%%%%%%%%%%%%%%%%%%%%%
% IMPORT SOLUTION Ux %
%%%%%%%%%%%%%%%%%%%%%%

tic;

rawData3 = importdata(file_ffSolUx);

solStrVect = rawData3.textdata(6:(end-1));
temp = strjoin(solStrVect);
ffSolUx = str2num(temp);

t3 = toc;

disp(['SOL. Ux - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t3),' [s]']);

%%%%%%%%%%%%%%%%%%%%%%
% IMPORT SOLUTION Uy %
%%%%%%%%%%%%%%%%%%%%%%

tic;

rawData4 = importdata(file_ffSolUy);

solStrVect = rawData4.textdata(6:(end-1));
temp = strjoin(solStrVect);
ffSolUy = str2num(temp);

t4 = toc;

disp(['SOL. Uy - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t4),' [s]']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT STIFFNESS MATRIX FOR PRESSURE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

rawData5 = importdata(file_Ap);

rows    = rawData5.data(2:end,1);
columns = rawData5.data(2:end,2);
values  = rawData5.data(2:end,3);

Ap = sparse(rows,columns,values);

t5 = toc;

disp(['STF. P  - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t5),' [s]']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT MASS MATRIX FOR PRESSURE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

rawData6 = importdata(file_Mp);

rows    = rawData6.data(2:end,1);
columns = rawData6.data(2:end,2);
values  = rawData6.data(2:end,3);

Mp = sparse(rows,columns,values);

t6 = toc;

disp(['MSS. P  - FINISHED EXTRACTING THE MESH in TIME t = ',num2str(t6),' [s]']);

end


