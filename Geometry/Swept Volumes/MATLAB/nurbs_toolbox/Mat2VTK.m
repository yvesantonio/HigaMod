function Mat2VTK(filename,matrix,format)
% Writes a 3D matrix as a *.VTK file as input for Paraview.
% Coded by Manuel A. Diaz, NHRI 08.21.2016

% Following the example VTK file:
% # vtk DataFile Version 2.0
% Volume example
% ASCII
% DATASET STRUCTURED_POINTS
% DIMENSIONS 3 4 6
% ASPECT_RATIO 1 1 1
% ORIGIN 0 0 0
% POINT_DATA 72
% SCALARS volume_scalars char 1
% LOOKUP_TABLE default
% 0 0 0 0 0 0 0 0 0 0 0 0
% 0 5 10 15 20 25 25 20 15 10 5 0
% 0 10 20 30 40 50 50 40 30 20 10 0
% 0 10 20 30 40 50 50 40 30 20 10 0
% 0 5 10 15 20 25 25 20 15 10 5 0
% 0 0 0 0 0 0 0 0 0 0 0 0

% Here the example is extended to n-dimensional matrices. e.g:
%
% matrix(:,:,1) = [0 0 0 0 0 0 0 0 0 0 0 0;
% 0 5 10 15 20 25 25 20 15 10 5 0;
% 0 10 20 30 40 50 50 40 30 20 10 0;
% 0 10 20 30 40 50 50 40 30 20 10 0;
% 0 5 10 15 20 25 25 20 15 10 5 0;
% 0 0 0 0 0 0 0 0 0 0 0 0];
% matrix(:,:,2) = [0 0 0 0 0 0 0 0 0 0 0 0;
% 0 5 10 15 20 25 25 20 15 10 5 0;
% 0 10 20 30 40 50 50 40 30 20 10 0;
% 0 10 20 30 40 50 50 40 30 20 10 0;
% 0 5 10 15 20 25 25 20 15 10 5 0;
% 0 0 0 0 0 0 0 0 0 0 0 0];
%
% Usage: Mat2VTK('example.vtk',matrix,'binary');

% Get the matrix dimensions.
[Nx,Ny,Nz] = size(matrix);

% Open the file.
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing.');
end

switch format
    case 'ascii'
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,'Volume example\n');
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET STRUCTURED_POINTS\n');
        fprintf(fid,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
        fprintf(fid,'ASPECT_RATIO %d %d %d\n',1,1,1);
        fprintf(fid,'ORIGIN %d %d %d\n',0,0,0);
        fprintf(fid,'POINT_DATA %d\n',Nx*Ny*Nz);
        fprintf(fid,'SCALARS Pressure int 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fwrite(fid, num2str(matrix(:)'));
    case 'binary'
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,'Volume example\n');
        fprintf(fid,'BINARY\n');
        fprintf(fid,'DATASET STRUCTURED_POINTS\n');
        fprintf(fid,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
        fprintf(fid,'ASPECT_RATIO %d %d %d\n',1,1,1);
        fprintf(fid,'ORIGIN %d %d %d\n',0,0,0);
        fprintf(fid,'POINT_DATA %d\n',Nx*Ny*Nz);
        fprintf(fid,'SCALARS Pressure float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fwrite(fid, matrix(:),'float','ieee-be');
    otherwise
        error('wrong input dummy :P');
end

% Close the file.
fclose(fid);
