clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%%
[v,f,n, name]=stlRead('grid sample.stl');
mystlPlot(v,f,0.3);hold on

x=1:0.5:19;
y=sqrt(x).*sin(x*0.5)+10;
z=zeros(numel(x),1);
plot3(x,y,z,'b')
s=[x' y' z];
%%
n=numel(s)/3;
list=zeros(n,1);
for i=1:n
    [xx yy zz index]=nearest(s(i,1),s(i,2),s(i,2),v(:,1),v(:,2),v(:,3));
    if norm([xx yy zz]-s(i,:))<1.e-5
        list(i)=index;
    else
        % find face which baricenter minimize distance point and baricenter
        j=nearest_face(s(i,:),f,v);
        % cancel j-th face and insert new vertex
        face=f(j,:);
        p=v(f(j,:),:);
        xgp=sum(p(:,1))/3;
        ygp=sum(p(:,2))/3;
        zgp=sum(p(:,3))/3;
        v(end+1,:)=[xgp ygp zgp];
        nf=zeros(numel(f)/3-1,3);  nf(1:j-1,:)=f(1:j-1,:); nf(j:end,:)=f(j+1:end,:);
        f=nf;
        f(end+1,:)=[face(1) face(2) numel(v)/3];
        f(end+1,:)=[face(2) face(3) numel(v)/3];
        f(end+1,:)=[face(3) face(1) numel(v)/3];
        list(i)=numel(v)/3;
    end      
end
linead=v(list,:);
mystlPlot(v,f,0.3);hold on
plot3(x,y,z,'b')
plot3(linead(:,1),linead(:,2),linead(:,3),'--r')
stlWrite('grid adapted.stl', f, v)
err=norm(linead-s)