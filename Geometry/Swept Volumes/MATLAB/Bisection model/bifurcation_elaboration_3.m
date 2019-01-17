%% define start and finish surface
clear all
close all
clc
load('data3')
s_start=surf_point1{1};
xbs=sum(s_start(:,1))/numel(s_start(:,1));ybs=sum(s_start(:,2))/numel(s_start(:,2));zbs=sum(s_start(:,3))/numel(s_start(:,3));
s_end=surf_point1{2};
xbe=sum(s_end(:,1))/numel(s_end(:,1));ybe=sum(s_end(:,2))/numel(s_end(:,2));zbe=sum(s_end(:,3))/numel(s_end(:,3));
s2=surf_point2{1};
xb=sum(s2(:,1))/numel(s2(:,1));yb=sum(s2(:,2))/numel(s2(:,2));zb=sum(s2(:,3))/numel(s2(:,3));
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
mystlPlot(vertices,faces,0.3);hold on
plot3(s2(:,1),s2(:,2),s2(:,3))
%% inserte meshes
% % n=numel(s2)/3;
% % list=zeros(n,1);
% % err=inf; errold=0;
% % for i=1:n
% %     [xx yy zz index]=nearest(s2(i,1),s2(i,2),s2(i,2),vertices(:,1),vertices(:,2),vertices(:,3));
% %     if norm([xx yy zz]-s2(i,:))<1
% %         list(i)=index;
% %     else
% %         % find face which baricenter minimize distance point and baricenter
% %         j=nearest_face(s2(i,:),faces,vertices);
% %         % cancel j-th face and insert new vertex
% %         face=faces(j,:);
% %         p=vertices(faces(j,:),:);
% %         xgp=sum(p(:,1))/3;
% %         ygp=sum(p(:,2))/3;
% %         zgp=sum(p(:,3))/3;
% %         vertices(end+1,:)=[xgp ygp zgp];
% %         nf=zeros(numel(faces)/3-1,3);  nf(1:j-1,:)=faces(1:j-1,:); nf(j:end,:)=faces(j+1:end,:);
% %         faces=nf;
% %         faces(end+1,:)=[face(1) face(2) numel(vertices)/3];
% %         faces(end+1,:)=[face(2) face(3) numel(vertices)/3];
% %         faces(end+1,:)=[face(3) face(1) numel(vertices)/3];
% %         list(i)=numel(vertices)/3;
% %     end      
% % end
% % linead=vertices(list,:);
% % mystlPlot(vertices,faces,0.3);hold on
% % plot3(s2(:,1),s2(:,2),s2(:,3),'b')
% % plot3(linead(:,1),linead(:,2),linead(:,3),'--r')
% % stlWrite('bifurcation grid adapted.stl', faces, vertices)
% % errold=err;
% % err=norm(linead-s2)
%% add pti curva sul stl
[x y z a]=nearest(s2(1,1),s2(1,2),s2(1,3),vertices(:,1),vertices(:,2),vertices(:,3));
p=0;
p(1)=a;
for i=2:numel(s2)/3
    [x y z a]=nearest(s2(i,1),s2(i,2),s2(i,3),vertices(:,1),vertices(:,2),vertices(:,3));
    if a~=p(end)
        p(end+1)=a;
    end
end
p=p';
xb=sum(s2(:,1))/numel(s2(:,1));yb=sum(s2(:,2))/numel(s2(:,2));zb=sum(s2(:,3))/numel(s2(:,3));
vertices(end+1,:)=[xb yb zb];
for j=1:numel(p)-1
    faces(end+1,:)=[p(j) p(j+1) numel(vertices)/3];    
    %normals(end+1,:)=mynormals(vertices(chained_edge{i}(j),:),vertices(chained_edge{i}(j+1),:),[xb yb zb]);    
end
    mystlPlot(vertices,faces,0.3);hold on
stlWrite('C:\Users\Leo\Documents\MATLAB\Bisection model\bifurcation_mod.stl',faces,vertices)
g=[xb yb zb]+TgC2(150,:);
hold on; plot3(g(1),g(2),g(3),'rd')
%% elliminare la terza edge ?????
 faces=edge_del(faces,vertices,[g]);
 mystlPlot(vertices,faces,0.3);hold on
stlWrite('C:\Users\Leo\Documents\MATLAB\Bisection model\bifurcation_mod.stl',faces,vertices)
hold on; plot3(g(1),g(2),g(3),'rd')
 
 %%
% f=faces;
% v=vertices;
% 
% i=1; n_faces=numel(f)/3;
% while  i<=n_faces
%     e=[f(i,1) f(i,2);
%        f(i,2) f(i,3); 
%        f(i,3) f(i,1)];
%    for j=1:3
%        u=1;
%        clearvars mult
%        n=numel(f)/3;
%        k=1;
%        while k<=n
%         if ( e(j,1)==f(k,1) | e(j,1)==f(k,2) | e(j,1)==f(k,3) ) & ...
%            ( e(j,2)==f(k,1) | e(j,2)==f(k,2) | e(j,2)==f(k,3) )
%                multi(u)=k;
%                u=u+1;
%         end
%         k=k+1;
%        end
%         if u>3
%             g1=sum( v(f(multi(1),:),:) )/3; d1=norm(g1-g); g2=sum( v(f(multi(2),:),:) )/3; d2=norm(g2-g); g3=sum( v(f(multi(3),:),:) )/3; d3=norm(g3-g);
%             fn=zeros(numel(f)/3-1,3);
%             if max([d1,d2,d3])==d3
%                 d=multi(3); fn(1:d-1,:)=f(1:d-1,:); 
%                 fn(d:end,:)=f(d+1:end,:); 
%                 f=fn;
%             end
%             if max([d1,d2,d3])==d2 
%                 d=multi(2); fn(1:d-1,:)=f(1:d-1,:); fn(d:end,:)=f(d+1:end,:); f=fn;
%             end
%             if max([d1,d2,d3])==d1  
%                 d=multi(1); fn(1:d-1,:)=f(1:d-1,:); fn(d:end,:)=f(d+1:end,:); f=fn;
%             end
%         end
%        end
%     n_faces=numel(f)/3;
%     i=i+1;
% end
%% dress edge
list=naked_edge(faces);
clearvars chained_edge
l=1;
while list~=-1
    [chained_edge{l}(:) list]=chain(list);
    l=l+1;
end
%% plot naked edge chained togheter
mystlPlot(vertices,faces,0.3);hold on
for i=1:numel(chained_edge)
plot3(vertices(chained_edge{i},1),vertices(chained_edge{i},2),vertices(chained_edge{i},3),'Color',[1 0 i/numel(chained_edge)],'LineWidth',3);hold on
end
% usare i baricentri per raccordare le centerline
%% create new faces from baricenter of linked edge 

xb=sum(vertices(chained_edge{1},1))/numel(vertices(chained_edge{1},1));
yb=sum(vertices(chained_edge{1},2))/numel(vertices(chained_edge{1},2));
zb=sum(vertices(chained_edge{1},3))/numel(vertices(chained_edge{1},3));
vertices(end+1,:)=[xb yb zb]-TgC2(150,:);
for j=1:numel(chained_edge{i})-1
    faces(end+1,:)=[chained_edge{i}(j:j+1) numel(vertices)/3];    
    %normals(end+1,:)=mynormals(vertices(chained_edge{i}(j),:),vertices(chained_edge{i}(j+1),:),[xb yb zb]);    
end

%% plot the dressed edge
mystlPlot(vertices,faces,0.3);hold on
plot3(vertices(end-1:end,1),vertices(end-1:end,2),vertices(end-1:end,3),'rd','MarkerSize',3);hold on
stlWrite('C:\Users\Leo\Documents\MATLAB\Bisection model\bifurcation_mod.stl',faces,vertices)
%% no problem mult has size 3 x 2 x ....  ; res=0
multi=edge_multeplicity(faces);
res=ismember(0,multi)