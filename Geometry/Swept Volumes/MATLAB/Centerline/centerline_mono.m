%% load data .dat
A=importdata('C:\Users\Leo\Desktop\cent_mono.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\model.stl');
% figure
% stlPlot(vertices,faces);hold on
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
%% get x y z
mystlPlot(vertices, faces, 0.1);
plot3(x,y,z,'LineWidth',3);grid on; xlabel('x'); ylabel('y')
