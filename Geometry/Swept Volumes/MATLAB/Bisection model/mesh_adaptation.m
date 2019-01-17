clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%%
load('data_3.mat')
[vertices,faces,nor, name]=stlRead('C:\Users\Leo\Desktop\bifurcation.stl');
mystlPlot(vertices,faces,0.3);hold on

s2=surf_point2{1};
plot3(s2(:,1),s2(:,2),s2(:,3),'r')
%%
n=numel(s2)/3;
list=zeros(n,1);
err=inf; errold=0;
while norm(err-errold)>1
for i=1:n
    [xx yy zz index]=nearest(s2(i,1),s2(i,2),s2(i,2),vertices(:,1),vertices(:,2),vertices(:,3));
    if norm([xx yy zz]-s2(i,:))<0.2
        list(i)=index;
    else
        % find face which baricenter minimize distance point and baricenter
        j=nearest_face(s2(i,:),faces,vertices);
        % cancel j-th face and insert new vertex
        face=faces(j,:);
        p=vertices(faces(j,:),:);
        xgp=sum(p(:,1))/3;
        ygp=sum(p(:,2))/3;
        zgp=sum(p(:,3))/3;
        vertices(end+1,:)=[xgp ygp zgp];
        nf=zeros(numel(faces)/3-1,3);  nf(1:j-1,:)=faces(1:j-1,:); nf(j:end,:)=faces(j+1:end,:);
        faces=nf;
        faces(end+1,:)=[face(1) face(2) numel(vertices)/3];
        faces(end+1,:)=[face(2) face(3) numel(vertices)/3];
        faces(end+1,:)=[face(3) face(1) numel(vertices)/3];
        list(i)=numel(vertices)/3;
    end      
end
linead=vertices(list,:);
mystlPlot(vertices,faces,0.3);hold on
plot3(s2(:,1),s2(:,2),s2(:,3),'b')
plot3(linead(:,1),linead(:,2),linead(:,3),'--r')
stlWrite('bifurcation grid adapted.stl', faces, vertices)
errold=err;
err=norm(linead-s2)
end