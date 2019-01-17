% %% sorting of the senterline to get separated curve
% mystlPlot(vertices,faces,0.3);hold on
% plot3(x,y,z)
% mystlPlot(vertices,faces,0.3);hold on
% for i=1:numel(x)
% plot3(x(i),y(i),z(i),'d','Color',[1-i/numel(x) 0 i/numel(x)])
% pause (0.05)
% end
clear all
close all
clc

load('data2.mat')
%% separete the 2 centerline
clearvars C1 C2 C3

C1(1,:)=[x(1) y(1) z(1)];
radius_C1(:,1)=radius_max(1);
jump=norm([x(1) y(1) z(1)]-[x(2) y(2) z(2)]);
k=2;
while jump<radius_max(k-1)
    C1(k,:)=[x(k) y(k) z(k)];
    radius_C1(k,1)=radius_max(k);
    jump=norm([x(k) y(k) z(k)]-[x(k+1) y(k+1) z(k+1)]);
    k=k+1;
end

not_in=1;
k=k+1;
i=2;
C2(1,:)=[x(k) y(k) z(k)];
radius_C2(1,1)=radius_max(k);
not_in=not_contained([x(k+1) y(k+1) z(k+1)],C1);
while not_in
    C2(i,:)=[x(k+1) y(k+1) z(k+1)];
    radius_C2(i,1)=radius_max(k+1);
    i=i+1;
    k=k+1;
    not_in=not_contained([x(k+1) y(k+1) z(k+1)],[C1]);
end
 
mystlPlot(vertices,faces,0.3);hold on
% plot3(x,y,z,'b')
plot3(C1(:,1),C1(:,2),C1(:,3),'--g',C2(:,1),C2(:,2),C2(:,3),'--r')
plot3(C1(1,1),C1(1,2),C1(1,3),'gd',C2(1,1),C2(1,2),C2(1,3),'rd')
%% 
mystlPlot(vertices, faces, 0.3);
%plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y'); hold on
[Cpt C]=nrb_approx(C1(:,1),C1(:,2),C1(:,3),15,3,1.e-10,30,100);
Cpt=nrbeval(C,[linspace(0,1,numel(Cpt)/3)]);
C1(:,1)=Cpt(1,:);
C1(:,2)=Cpt(2,:);
C1(:,3)=Cpt(3,:);
clearvars Cpt

plot3(C1(:,1),C1(:,2),C1(:,3),'g');grid on; xlabel('x'); ylabel('y')
plot3(C1(1,1),C1(1,2),C1(1,3),'rd');grid on; xlabel('x'); ylabel('y');
plot3(C1(end,1),C1(end,2),C1(end,3),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

%%

%% plot vector tg nm bn

% % DEBUG: PLOT THE VECTOR TANGENT TO THE CENTERLINE IN FIVE DIFFERENT
% % LOCATIONS
% step=1;
% figure;
% hold on;
% for i = 1:step:numel(x)
%     pInit=[x(i) y(i) z(i)];
%     pFinal=[x(i)+Tg(i,1) y(i)+Tg(i,2) z(i)+Tg(i,3)];      
%     vectarrow(pInit,pFinal);
%     hold on;
% end
% for i = 1:step:numel(x)
%     pInit=[x(i) y(i) z(i)];
%     pFinal=[x(i)+Bn(i,1) y(i)+Bn(i,2) z(i)+Bn(i,3)];      
%     vectarrow(pInit,pFinal);
%     
%     hold on;
% end
% for i = 1:step:numel(x)
%     pInit=[x(i) y(i) z(i)];
%     pFinal=[x(i)+Nm(i,1) y(i)+Nm(i,2) z(i)+Nm(i,3)];      
%     vectarrow(pInit,pFinal);
%     
%     hold on;
% end
% hold off;
% 
%% DA VEDERE!!! 
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
% plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
% plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')
% prima curva
k=1;
step=2;
[TgC1 NmC1 BnC1]=myfrenet(C1(:,1),C1(:,2),C1(:,3));
for i=[170 195]
    i
   plane = [
        C1(i,1),C1(i,2), C1(i,3), ...
        NmC1(i,:), ...
        BnC1(i,:)
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
        surf_point1{k}=[polyCellNow{min(1)}(:,1) polyCellNow{min(1)}(:,2) polyCellNow{min(1)}(:,3)]   ;
        x_plot=polyCellNow{min(1)}(:,1)';  y_plot=polyCellNow{min(1)}(:,2)';  z_plot=polyCellNow{min(1)}(:,3)';
        plot3(x_plot(1,:),y_plot(1,:),z_plot(1,:),'Color',[1-i/numel(x) 0 i/numel(x)],'LineWidth',0.1);grid on;hold on 
        k=k+1;
    end
pause(0.01)
end 
%% seconda curva
[Cpt C]=nrb_approx(C2(:,1),C2(:,2),C2(:,3),14,3,1.e-10,30,100);
Cpt=nrbeval(C,[linspace(0,1,numel(Cpt)/3)]);
C2(:,1)=Cpt(1,:);
C2(:,2)=Cpt(2,:);
C2(:,3)=Cpt(3,:);
clearvars Cpt
[TgC2 NmC2 BnC2]=myfrenet(C2(:,1),C2(:,2),C2(:,3));

plot3(C2(:,1),C2(:,2),C2(:,3),'r');grid on; xlabel('x'); ylabel('y')
k=1;
step=2;
for i=150
    i
   plane = [
        C2(i,1), C2(i,2), C2(i,3), ...
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
        surf_point2{k}=[polyCellNow{min(1)}(:,1) polyCellNow{min(1)}(:,2) polyCellNow{min(1)}(:,3)]   ;
        x_plot=polyCellNow{min(1)}(:,1)';  y_plot=polyCellNow{min(1)}(:,2)';  z_plot=polyCellNow{min(1)}(:,3)';
        plot3(x_plot(1,:),y_plot(1,:),z_plot(1,:),'Color',[ 0 1-i/numel(x) i/numel(x)],'LineWidth',0.1);grid on;hold on        
        k=k+1;
    end
pause(0.01)
end 
%%

    