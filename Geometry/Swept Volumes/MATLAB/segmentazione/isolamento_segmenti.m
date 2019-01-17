clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% load data .dat , need to run the vmtk to get the .dat file

[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\marcione.stl');


%%
mystlPlot(vertices,faces,0.3), xlabel('x'),ylabel('y'); hold on

 plane = [64 2.5 52, ...
          0 0 1, ...
          1 0 0];
polyCellNow = xsecmesh(plane, vertices, faces);
sec1=[polyCellNow{1}(:,1) polyCellNow{1}(:,2) polyCellNow{1}(:,3)];

%  plane = [0 1.25 0, ...
%           0 0 1, ...
%           1 0 0];
% polyCellNow = xsecmesh(plane, vertices, faces);
% sec2=[polyCellNow{1}(:,1) polyCellNow{1}(:,2) polyCellNow{1}(:,3)];

plot3(sec1(:,1),sec1(:,2),sec1(:,3),'rd')
% plot3(sec2(:,1),sec2(:,2),sec2(:,3),'bd')

%%
load('data_marcione_modificato.mat')
mystlPlot(vertices,faces,0.3), xlabel('x'),ylabel('y'); hold on

%%



[x y z a]=nearest(sec1(1,1),sec1(1,2),sec1(1,3),vertices(:,1),vertices(:,2),vertices(:,3));
p=0;
p(1)=a;
for i=2:numel(sec1)/3
    [x y z a]=nearest(sec1(i,1),sec1(i,2),sec1(i,3),vertices(:,1),vertices(:,2),vertices(:,3));
    if a~=p(end)
        p(end+1)=a;
    end
end
p=p';
plot3(vertices(p,1),vertices(p,2),vertices(p,3),'g','LineWidth',2)
for i=1:numel(p)-1
    limit(i,:)=[p(i) p(i+1)];
end
%% isolamento facce
f=faces;

mystlPlot(vertices,f,0.3)
ne=find_f_naked(f);
plot3(vertices(ne(1:2),1),vertices(ne(1:2),2),vertices(ne(1:2),3),'b','LineWidth',2)
nf1=[];
[nf1 f]=segmento(nf1,f,ne(1:2),ne(3),St,vertices);
mystlPlot(vertices,nf1,0.3)
pause(1)
disp('diocan');
%%
f=newf(f);
ne=find_f_naked(f);
plot3(vertices(ne(1:2),1),vertices(ne(1:2),2),vertices(ne(1:2),3),'b','LineWidth',2)
nf2=[];
[nf2 f]=segmento(nf2,f,ne(1:2),ne(3),St,vertices);
mystlPlot(vertices,nf2,0.3)
pause(1)

    
    