%% definition of tube, " noise added "
clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
% vmtkcenterlines -ifile C:/Users/Leo/Desktop/my_model.stl  -ofile C:/Users/Leo/Desktop/cent_my_model.dat 
%% load data .dat , need to run the vmtk to get the .dat file
A=importdata('C:\Users\Leo\Desktop\cent_my_model.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\my_model.stl');
stlPlot(vertices,faces,name)
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
%% dress edge after check that that exist
list=naked_edge(faces);
if list~=[-1 -1]

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
        vertices(end+1,:)=[xb yb zb];
        for j=1:numel(chained_edge{i})-1
            faces(end+1,:)=[chained_edge{i}(j:j+1) numel(vertices)/3];    
            normals(end+1,:)=mynormals(vertices(chained_edge{i}(j),:),vertices(chained_edge{i}(j+1),:),[xb yb zb]);    
        end
        end
        %% plot the dressed edge
        mystlPlot(vertices,faces,1);hold on
        plot3(vertices(end-1:end,1),vertices(end-1:end,2),vertices(end-1:end,3),'rd','MarkerSize',3);hold on
end
%%
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y'); hold on
[Cpt C]=nrb_approx(x,y,z,15,3,1.e-10,30,100);
Cpt=nrbeval(C,[linspace(0,1,numel(Cpt)/3)]);
x=Cpt(1,:);
y=Cpt(2,:);
z=Cpt(3,:);
clearvars Cpt

plot3(x,y,z,'g');grid on; xlabel('x'); ylabel('y')
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
step=1;
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

%% DA VEDERE!!! 
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

k=1;
step=2;
for i=[1:10:200]
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
        
        coordtg(k,:)=[x(i) y(i) z(i)];
        tgsrf(k,:)=Tg(i,:);
        k=k+1;
    end
pause(0.01)
end 
save('C:\Users\Leo\Documents\MATLAB\Volume\section_model_step2.mat','surf_point')
save('C:\Users\Leo\Documents\MATLAB\Volume\Tg_step2.mat','tgsrf')
save('C:\Users\Leo\Documents\MATLAB\Volume\coord_step2.mat','coordtg')



