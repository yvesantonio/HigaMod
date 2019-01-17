
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

plot3(sec1(:,1),sec1(:,2),sec1(:,3),'rd')
% plot3(sec2(:,1),sec2(:,2),sec2(:,3),'bd')


%% inserisci le mesh che servono


for j=1:numel(sec1)/3-1
    p1=sec1(j,:);
    p2=sec1(j+1,:);
i=1;
while i 
   %cerca la faccia che contenga il punto j-esimo e quelli j+1-esimo
   f=faces(i,:);

    v=[vertices(f(1),:); vertices(f(2),:); vertices(f(3),:)];
    index=[1 2; 2 3; 3 1];
    existp1=0;
    existp2=0;
    for k=1:3
        t=point_in_edge(p1,v(index(k,1),:),v(index(k,2),:));
        if t==1
            existp1=k;
        end
    end
    if existp1~=0
    for k=1:3
        t=point_in_edge(p2,v(index(k,1),:),v(index(k,2),:));
        if t==1
            existp2=k;
        end
    end
    end
    
if existp1~=0 & existp2~=0
    
% look for common vertex of edge
e1=[f(index(existp1,:))];
e2=[f(index(existp2,:))];

if ismember(e1(1),e2)
    cv=e1(1);  pnc1=e1(2);
else
   cv=e1(2);  pnc1=e1(1);
end

if ismember(e2(1),e1)
    cvt=e2(1);  pnc2=e2(2);
else
   cvt=e2(2);  pnc2=e2(1);
end


        break
    end

    i=i+1;   
end

[vertices ve1]=addvertices(vertices,p1);
[vertices ve2]=addvertices(vertices,p2);
faces(i,:)=[ve1 ve2 cv];
faces(end+1,:)=[ve1 ve2 pnc1];
faces(end+1,:)=[ve2 pnc2 pnc1];
end
%%
for j=1:numel(  index_face  )-1
    faces(index_face(j),:)=[ve1(j) ve2(j) cv(j)];
    faces(end+1,:)=[ve1(j) ve2(j) pnc1(j)];
    faces(end+1,:)=[ve2(j) pnc2(j) pnc1(j)];
end
mystlPlot(vertices,faces,1), xlabel('x'),ylabel('y'); hold on

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
% [xn yn zn]=nearest(xb,yb,zb,x,y,z);
% if [xn yn zn]==[x(1) y(1) z(1)]
%     x(2:end+1)=x; x(1)=xb;
%     y(2:end+1)=y; y(1)=yb;
%     z(2:end+1)=z; z(1)=zb;
% else
%     x(end+1)=xb;
%     y(end+1)=yb;
%     z(end+1)=zb;    
% end
vertices(end+1,:)=[xb yb zb];
for j=1:numel(chained_edge{i})-1
    faces(end+1,:)=[chained_edge{i}(j:j+1) numel(vertices)/3];    
    normals(end+1,:)=mynormals(vertices(chained_edge{i}(j),:),vertices(chained_edge{i}(j+1),:),[xb yb zb]);    
end
end

mystlPlot(vertices,faces,0.1), xlabel('x'),ylabel('y'); hold on





