clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
A=importdata('C:\Users\Leo\Desktop\cent_bifurcation.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\bifurcation.stl');
stlPlot(vertices,faces,name); xlabel('x');ylabel('y');zlabel('z')
%%
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
radius_max=A.data(:,4);
%% dress edge
list=naked_edge(faces);
l=1;
while list~=-1
    [chained_edge{l}(:) list]=chain(list);
    l=l+1;
end
%% plot naked edge chained togheter
stlPlot(vertices,faces,name);hold on
for i=1:numel(chained_edge)
plot3(vertices(chained_edge{i},1),vertices(chained_edge{i},2),vertices(chained_edge{i},3),'Color',[1 0 i/numel(chained_edge)],'LineWidth',3);hold on
end
% usare i baricentri per raccordare le centerline
%% create new faces from baricenter of linked edge 
for i=1:numel(chained_edge)
xb=sum(vertices(chained_edge{i},1))/numel(vertices(chained_edge{i},1));
yb=sum(vertices(chained_edge{i},2))/numel(vertices(chained_edge{i},2));
zb=sum(vertices(chained_edge{i},3))/numel(vertices(chained_edge{i},3));
vertices(end+1,:)=[xb yb zb];
for j=1:numel(chained_edge{i})-1
    faces(end+1,:)=[chained_edge{i}(j:j+1) numel(vertices)/3];    
    normals(end+1,:)=mynormals(vertices(chained_edge{i}(j),:),vertices(chained_edge{i}(j+1),:),[xb yb zb]);    
end
end
%% plot the dressed edge
mystlPlot(vertices,faces,0.3);hold on
plot3(vertices(end-1:end,1),vertices(end-1:end,2),vertices(end-1:end,3),'rd','MarkerSize',3);hold on
stlWrite('bifurcation_closed.stl',faces,vertices)


