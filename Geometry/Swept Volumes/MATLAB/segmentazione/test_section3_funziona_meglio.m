clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% load data .dat , need to run the vmtk to get the .dat file

[v,f,normals,name] = stlRead('C:\Users\Leo\Desktop\bifurcation.stl');
mystlPlot(v,f,1), xlabel('x'),ylabel('y'); hold on

 plane = [64 2.5 52, ...
          0 0 1, ...
          1 0 0];
polyCellNow = xsecmesh(plane, v, f);
sec=[polyCellNow{1}(:,1) polyCellNow{1}(:,2) polyCellNow{1}(:,3)];

%  plane = [0 1.25 0, ...
%           0 0 1, ...
%           1 0 0];
% polyCellNow = xsecmesh(plane, vertices, faces);
% sec2=[polyCellNow{1}(:,1) polyCellNow{1}(:,2) polyCellNow{1}(:,3)];


 plot3(sec(17,1),sec(17,2),sec(17,3),'bd')
%%
index=[1 2;2 3; 3 1];
for n=1:numel(sec)/3
    p=sec(n,:);
    dist(n,1)=inf;
for i=1:numel(f)/3
    fi=f(i,:);
    for ii=1:3
        v1=v(fi(index(ii,1)),:); v2=v(fi(index(ii,2)),:);
        d=abs(v2-v1);
        
%         if prod(abs(v2-p)<d & abs(v1-p)<d)
            dist_new=norm(cross(p-v1,v2-v1))/norm(v2-v1);
            if dist_new<dist(n,1)
                dist(n,1)=dist_new;
                edge(n,:)=fi(index(ii,:));
                plot3([v1(1) v2(1)],[v1(2) v2(2)] ,[v1(3) v2(3)] ,'r')
%             end
        
        end
    end
end
end
faccia=zeros(size(edge));
for i=1:numel(sec)/3
    k=1;
    for ii=1:numel(f)/3
        if ismember(edge(i,:),f(ii,:))
            faccia(i,k)=ii;
            k=k+1;
            if k==3
                break
            end
        end
    end
end
%% inserimento mesh
for i=1:numel(sec)/3-1
    
    if ismember(faccia(i,1),faccia(i+1,:))
        cf=faccia(i,1);
    else
        cf=faccia(i,2);
    end
    
    if ismember(edge(i,1),edge(i+1,:))
        cv=edge(i,1); ncv1=edge(i,2);
    else
        cv=edge(i,2); ncv1=edge(i,1);
    end
    
    if ismember(edge(i+1,1),edge(i,:))
        ncv2=edge(i+1,2);
    else
        ncv2=edge(i+1,1);
    end
    [v, vp1]=addvertices(v,sec(i,:));
    [v, vp2]=addvertices(v,sec(i+1,:));
    f(cf,:)=[cv vp1 vp2];
    f(end+1,:)=[vp1 vp2 ncv1];
    f(end+1,:)=[ncv1 ncv2 vp2];
    
    
    
    
    
    
    
end

mystlPlot(v, f, 1);hold on 
plot3(sec(:,1),sec(:,2),sec(:,3),'--g')
multi=edge_multeplicity(f);  



