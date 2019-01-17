%%
clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
load('data4');
%%
mystlPlot(vertices, faces, 0.3); hold on;grid on;
% plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
% plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
% plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')
% prima curva
%surf_point{1}=s_start;
k=1;
for i=[100:1:230]
    i
   plane = [
            C1(i,1),C1(i,2), C1(i,3), ...
            NmC1(i,:), ...
            BnC1(i,:)] ;
    % Get the plane/mesh cross-section.
    polyCellNow = xsecmesh(plane, vertices, faces);
    % Store the results.
    if numel(polyCellNow)~=0
            min=[0 inf];
            for j=1:numel(polyCellNow)
                xm=sum(polyCellNow{j}(:,1))/numel(polyCellNow{j}(:,1));
                ym=sum(polyCellNow{j}(:,2))/numel(polyCellNow{j}(:,2));
                zm=sum(polyCellNow{j}(:,3))/numel(polyCellNow{j}(:,3));
                if norm([xm-C1(i,1) ym-C1(i,2) zm-C1(i,3)])<min(2)
                    min=[j norm([xm-C1(i,1) ym-C1(i,2) zm-C1(i,3)])];
                end
            end
        surf_point_main{k}=[polyCellNow{min(1)}(:,1) polyCellNow{min(1)}(:,2) polyCellNow{min(1)}(:,3)]   ;
        coordtg_main(k,:)=C1(i,:);
        Tg_section_main(k,:)=TgC1(i,:);
        plot3(surf_point_main{k}(:,1),surf_point_main{k}(:,2),surf_point_main{k}(:,3),'b')
        k=k+1;
    end
pause(0.01)
end 
save('data_surf_main','surf_point_main','coordtg_main','Tg_section_main')
%surf_point{end+1}=s_end;

%%
[vertices,faces,normals,name] = stlRead('bifurcation_closed.stl');
clearvars coordtg Tg_section surf_point_sec
k=1;
step=2;
[TgC2 NmC2 BnC2]=myfrenet(C2(:,1),C2(:,2),C2(:,3));
for i=[20:1:150]
    i
   plane = [
        C2(i,1),C2(i,2), C2(i,3), ...
        NmC2(i,:), ...
        BnC2(i,:)];
    % Get the plane/mesh cross-section.
    polyCellNow = xsecmesh(plane, vertices, faces);
    % Store the results.
    if numel(polyCellNow)~=0
            min=[0 inf];
            for j=1:numel(polyCellNow)
                xm=sum(polyCellNow{j}(:,1))/numel(polyCellNow{j}(:,1));
                ym=sum(polyCellNow{j}(:,2))/numel(polyCellNow{j}(:,2));
                zm=sum(polyCellNow{j}(:,3))/numel(polyCellNow{j}(:,3));
                if norm([xm-C2(i,1) ym-C2(i,2) zm-C2(i,3)])<min(2)
                    min=[j norm([xm-C2(i,1) ym-C2(i,2) zm-C2(i,3)])];
                end
            end
        surf_point_sec{k}=[polyCellNow{min(1)}(:,1) polyCellNow{min(1)}(:,2) polyCellNow{min(1)}(:,3)]   ;
        coordtg_sec(k,:)=C2(i,:);
        Tg_section_sec(k,:)=TgC2(i,:);
        plot3(surf_point_sec{k}(:,1),surf_point_sec{k}(:,2),surf_point_sec{k}(:,3))
        k=k+1;
    end
pause(0.01)
end 
save('data_surf_sec','surf_point_sec','coordtg_sec','Tg_section_sec')


