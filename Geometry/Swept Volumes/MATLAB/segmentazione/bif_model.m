clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
load('bif_2_completa.mat')
bif=St.bif;
c1=St.next1.data(St.fine:St.next1.start,:);
c2=St.next2.data(St.fine:St.next2.start,:);
mystlPlot(v, bif, 0.3);hold on
plot3(c1(:,1),c1(:,2),c1(:,3),'b','LineWidth',2);plot3(c1(1,1),c1(1,2),c1(1,3),'bd')
plot3(c2(:,1),c2(:,2),c2(:,3),'r','LineWidth',2);plot3(c2(1,1),c2(1,2),c2(1,3),'rd')
%%
[faces v]=closefaces(bif,v);
k=1;
for i=[10 20 30 45]
    j=i;
    [ srf_pti1{k} j ct t]=section_mesh(c1, v, faces, j);
    plot3(srf_pti1{k}(:,1),srf_pti1{k}(:,2),srf_pti1{k}(:,3),'b','LineWidth',2)
    k=k+1;
end
k=1;
for i=[10 20 30 45]
    j=i;
    [ srf_pti2{k} j ct t]=section_mesh(c2, v, faces, j);
    plot3(srf_pti2{k}(:,1),srf_pti2{k}(:,2),srf_pti2{k}(:,3),'r','LineWidth',2)
    k=k+1;
end