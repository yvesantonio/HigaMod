%% definition of tube, " noise added "
clear all
close all
clc
%% load data .dat , need to run the vmtk to get the .dat file
A=importdata('C:\Users\Leo\Desktop\cent_my_model.dat'); % dati relativi alla centerline
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\my_model.stl');
stlPlot(vertices,faces,name)
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
%% 
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y'); hold on
[Cpt C]=nrb_approx(x,y,z,30,3,1.e-10,30,100); % nurbs fitting for continuity
Cpt=nrbeval(C,[linspace(0,1,2*numel(Cpt))]); % calcola tanti punti equidistanti per una migliore my_frenet
x=Cpt(1,:);
y=Cpt(2,:);
z=Cpt(3,:);
clearvars Cpt

plot3(x,y,z,'g');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

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

%% non necessiario fare run, ho già i dati pronti
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

k=1;
step=6;
for i=[1:step:numel(x)]
    i
   plane = [
        x(i), y(i), z(i), ...
        Nm(i,:), ...
        Bn(i,:)
        ];
    % Get the plane/mesh cross-section.
    polyCellNow = xsecmesh(plane, vertices, faces);
    % Store the results.
    if numel(polyCellNow)~=0
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
        plot3(x_plot(1,:),y_plot(1,:),z_plot(1,:),'Color',[1-i/numel(x) 0 i/numel(x)],'LineWidth',0.1);grid on;hold on 
        k=k+1;
    end
pause(1)
end 
save('C:\Users\Leo\Documents\MATLAB\Curve fitting\model_to_section_yves.mat','surf_point')

% vedi figura " model_to_section_yves "

%% design the coons surface related to each section
% [vertices,faces,normals,name] = stlRead('my_model.stl'); 
hold on
stlPlot(vertices,faces,0.3)
load('model_to_section_yves.mat')
for i=[1 10:5:numel(surf_point) numel(surf_point)]
    coons_srf(i)=get_coons(surf_point{i},6,2);
   
    nrbplot(coons_srf(i),[10 10])
end