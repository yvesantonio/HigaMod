%% create the junction
clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
Data=importdata('C:\Users\Leo\Desktop\cent_bifurcation.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\bifurcation.stl');
load('surf1.mat')
load('surf2.mat')
%%
s1=surf_point1{1};s2=surf_point1{2};s3=surf_point2{1};
bar(1,:)=[sum(s1(:,1))/numel(s1(:,1)) sum(s1(:,2))/numel(s1(:,1)) sum(s1(:,3))/numel(s1(:,1))];
bar(2,:)=[sum(s2(:,1))/numel(s2(:,1)) sum(s2(:,2))/numel(s2(:,1)) sum(s2(:,3))/numel(s2(:,1))];
bar(3,:)=[sum(s3(:,1))/numel(s3(:,1)) sum(s3(:,2))/numel(s3(:,1)) sum(s3(:,3))/numel(s3(:,1))];
bar(4,:)=[sum(bar(1:3,1))/3 sum(bar(1:3,2))/3 sum(bar(1:3,3))/3 ];
mystlPlot(vertices, faces, 0.3); hold on;
plot3(s1(:,1),s1(:,2),s1(:,3),'b',...
      s2(:,1),s2(:,2),s2(:,3),'b',...
      s3(:,1),s3(:,2),s3(:,3),'b',...
      bar(4,1),bar(4,2),bar(4,3),'bd');grid on; xlabel('x'); ylabel('y')
v(1,:)=[bar(4,:)]-[bar(1,:)]; v(2,:)=[bar(4,:)]-[bar(2,:)]; v(3,:)=[bar(4,:)]-[bar(3,:)];
vectarrow(bar(1,:),bar(1,:)+v(1,:)); hold on;
vectarrow(bar(2,:),bar(4,:)); hold on;
vectarrow(bar(3,:),bar(4,:)); hold on;
%% definie piano e vettore ad esso perpendicolare
perp=mygaus(v,[0;0;0],[1;1;1],1.e-2,100); perp=perp/norm(perp);
vectarrow(bar(4,:),bar(4,:)+perp'); hold on;