%% definition of tube, " noise added "
clear all
close all
clc
%% load data .dat , need to run the vmtk to get the .dat file
A=importdata('C:\Users\Leo\Desktop\cent_model.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\model.stl');
stlPlot(vertices,faces,name)
%% get x y z 
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);

%% filtraggio
b=ones(1,10)/10; 
xfil = filter(b,1,x);
yfil = filter(b,1,y);
zfil = filter(b,1,z);
plot3(x,y,z,'b',xfil(numel(b):end),yfil(numel(b):end),zfil(numel(b):end),'r');xlabel('x');ylabel('y')
                                                                       ;title('r filtered, b original')
x=xfil(numel(b):end);
y=yfil(numel(b):end);
z=zfil(numel(b):end);
mystlPlot(vertices, faces, 0.1);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

%%
[Tg Nm Bn]=myfrenet(x,y,z);

%%

mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd',x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')


k=1;
for i=1:2:180 %numel(x)
    i
   plane = [
        x(i), y(i), z(i), ...
        Bn(i,:), ...
        Nm(i,:)
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
plot3(x_plot(1,:),y_plot(1,:),z_plot(1,:),'Color',[0.5 0 i/numel(x)],'LineWidth',0.1);grid on;hold on 
k=k+1;
end










    