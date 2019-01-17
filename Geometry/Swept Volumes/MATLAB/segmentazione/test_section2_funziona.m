clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% load data .dat , need to run the vmtk to get the .dat file

[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\marcione.stl');
mystlPlot(vertices,faces,1), xlabel('x'),ylabel('y'); hold on

 plane = [64 2.75 52, ...
          0 0 1, ...
          1 0 0];
polyCellNow = xsecmesh(plane, vertices, faces);
sec1=[polyCellNow{1}(:,1) polyCellNow{1}(:,2) polyCellNow{1}(:,3)];

%  plane = [0 1.25 0, ...
%           0 0 1, ...
%           1 0 0];
% polyCellNow = xsecmesh(plane, vertices, faces);
% sec2=[polyCellNow{1}(:,1) polyCellNow{1}(:,2) polyCellNow{1}(:,3)];


% plot3(sec2(:,1),sec2(:,2),sec2(:,3),'bd')
%%
j=1;
nsec=[];
for i=1:numel(sec1)/3
    [x y z index d(i,1)]=mynearest(sec1(i,1),sec1(i,2),sec1(i,3),vertices(:,1),vertices(:,2),vertices(:,3));
    if d(i)<1.e-2
        is=isin(vertices(index,1:3),nsec(:,1:3));
        if ~is
        nsec(j,:)=[vertices(index,:) 1];
        j=j+1;
        end
       
    else
        nsec(j,:)=[sec1(i,:) 0];
        j=j+1;
    end
    
end
plot3(nsec(:,1),nsec(:,2),nsec(:,3),'g','LineWidth',2)
p=nsec;
clearvars index
%%
mystlPlot(vertices,faces,1)
for i=1:numel(p)/4-1
p1=p(i,:); p2=p(i+1,:);

    if p1(4)==0 & p2(4)==0
[faces vertices]=addfaces00(faces,vertices,p1(1:3),p2(1:3));

    end
    
    if p1(4)==1 & p2(4)==0
[faces vertices]=addfaces10(faces,vertices,p1(1:3),p2(1:3));

plot3(nsec(i:i+1,1),nsec(16:17,2),nsec(i:i+1,3),'go')
    end
    
    
    if p1(4)==0 & p2(4)==1
            [faces vertices]=addfaces10(faces,vertices,p2(1:3),p1(1:3));

plot3(nsec(i:i+1,1),nsec(16:17,2),nsec(i:i+1,3),'go')
     end
    
        
end
mystlPlot(vertices,faces,1)
    
    
  



