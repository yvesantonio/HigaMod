geometry = 'Slab';

% fileToRead1 = [pwd,'/Geometry/3D_',geometry,'/ff/ffMesh.mesh'];
% fileToRead2 = [pwd,'/Geometry/3D_',geometry,'/ff/ffSolution.sol'];
% fileToRead3 = [pwd,'/Geometry/3D_',geometry,'/ff/A.txt'];
% fileToRead4 = [pwd,'/Geometry/3D_',geometry,'/ff/M.txt'];

fileToRead1 = 'ffMesh.mesh';
fileToRead2 = 'ffSolution.sol';
fileToRead3 = 'A.txt';
fileToRead4 = 'M.txt';


[structMesh,ffSol,stiffMat,massMat] = my_FFimportfilemesh_3D(fileToRead1, ...   
                                                             fileToRead2, ...
                                                             fileToRead3, ...
                                                             fileToRead4);
                                                         
filename1 = 'ffSol.mat';
filename2 = 'stiffMat.mat';
filename3 = 'massMat.mat';
filename4 = 'structMesh.mat';

save(filename1,'ffSol');
save(filename2,'stiffMat');
save(filename3,'massMat');
save(filename4,'structMesh');