%% definition of tube, " noise added "
clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%%
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\pipe.stl');
sp1=linspace(0,10*pi,numel(vertices(:,1)))';
sp2=linspace(0,3*pi,numel(vertices(:,1)))';
vertices(:,3)=vertices(:,3)+randi([-1000 1000],numel(vertices(:,1)),1)*5.*sin(sp1)*1.e-5;                        
%vertices(:,2)=vertices(:,2)+randi([-1000 1000],numel(vertices(:,1)),1).*sin(sp1)*0.5*1.e-4;
%vertices(:,1)=vertices(:,1)+randi([-1000 1000],numel(vertices(:,1)),1).*sin(sp1)*0.5*1.e-4;
stlPlot(vertices,faces,name);
stlWrite('C:\Users\Leo\Desktop\pipe_mod.stl',faces,vertices)
%% load data .dat , need to run the vmtk to get the .dat file
A=importdata('C:\Users\Leo\Desktop\cent_pipe.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\pipe_mod.stl');
stlPlot(vertices,faces,name)
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
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
%% plot the dressed edge
mystlPlot(vertices,faces,0.3);hold on
plot3(vertices(end-1:end,1),vertices(end-1:end,2),vertices(end-1:end,3),'rd','MarkerSize',3);hold on
%%
[Cpt C]=nrb_approx(x,y,z,10,2,1.e-7,30,100);
x=Cpt(1,:);
y=Cpt(2,:);
z=Cpt(3,:);
clearvars Cpt
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

%% filtraggio
% b=ones(1,30)/30; 
% xfil = filter(b,1,x);
% yfil = filter(b,1,y);
% zfil = filter(b,1,z);
% figure
% plot3(x,y,z,'b',xfil(numel(b):end),yfil(numel(b):end),zfil(numel(b):end),'r');xlabel('x');ylabel('y');grid on
%                                                                             ;title('r filtered, b original')
% x=xfil(numel(b):end);
% y=yfil(numel(b):end);
% z=zfil(numel(b):end);

%%
[Tg Nm Bn]=myfrenet(x,y,z);
%% plot vector tg nm bn

% % DEBUG: PLOT THE VECTOR TANGENT TO THE CENTERLINE IN FIVE DIFFERENT
% % LOCATIONS
step=5;
figure;
hold on;
for i = 1:step:numel(x)
    pInit=[x(i) y(i) z(i)];
    pFinal=[x(i)+Tg(i,1) y(i)+Tg(i,2) z(i)+Tg(i,3)];      
    vectarrow(pInit,pFinal);
    hold on;
end
for i = 1:step:numel(x)
    pInit=[x(i) y(i) z(i)];
    pFinal=[x(i)+Bn(i,1) y(i)+Bn(i,2) z(i)+Bn(i,3)];      
    vectarrow(pInit,pFinal);
    
    hold on;
end
for i = 1:step:numel(x)
    pInit=[x(i) y(i) z(i)];
    pFinal=[x(i)+Nm(i,1) y(i)+Nm(i,2) z(i)+Nm(i,3)];      
    vectarrow(pInit,pFinal);
    
    hold on;
end
hold off;

%%
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

surf_point{1}=[vertices(chained_edge{2},1) vertices(chained_edge{2},2) vertices(chained_edge{2},3)];
x_plot=surf_point{1}(:,1);  y_plot=surf_point{1}(:,2);  z_plot=surf_point{1}(:,3);
plot3(x_plot,y_plot,z_plot,'Color',[1 0 0],'LineWidth',2);grid on;hold on 
k=2;
step=100;
for i=1:step:numel(x)
    i
   plane = [
        x(i), y(i), z(i), ...
        Nm(i,:), ...
        Bn(i,:)
        ];
    % Get the plane/mesh cross-section.
    polyCellNow = xsecmesh(plane, vertices, faces);
    % Store the results.
    min=[0 inf];
    for j=1:numel(polyCellNow)
        xm=sum(polyCellNow{j}(:,1))/numel(polyCellNow{j}(:,1));
        ym=sum(polyCellNow{j}(:,2))/numel(polyCellNow{j}(:,2));
        zm=sum(polyCellNow{j}(:,3))/numel(polyCellNow{j}(:,3));
        if norm([xm-x(i) ym-y(i) zm-z(i)])<min(2)
            min=[j norm([xm-x(i) ym-y(i) zm-z(i)])];
        end
    end
surf_point{k}=[polyCellNow{min(1)}(:,1) polyCellNow{min(1)}(:,2) polyCellNow{min(1)}(:,3)]   ;
x_plot=polyCellNow{min(1)}(:,1)';  y_plot=polyCellNow{min(1)}(:,2)';  z_plot=polyCellNow{min(1)}(:,3)';
plot3(x_plot(1,:),y_plot(1,:),z_plot(1,:),'Color',[1-i/numel(x) 0 i/numel(x)],'LineWidth',2);grid on;hold on 
k=k+1;

end 
surf_point{k}=[vertices(chained_edge{1},1) vertices(chained_edge{1},2) vertices(chained_edge{1},3)]   ;
x_plot=surf_point{k}(:,1);  y_plot=surf_point{k}(:,2);  z_plot=surf_point{k}(:,3);
plot3(x_plot,y_plot,z_plot,'Color',[0 0 1],'LineWidth',2);grid on;hold on ;title('red start, blue end')
