%% tubo troll
y=0:0.1:15;
x=zeros(151,1);
z=x;
%A=importdata('C:\Users\Leo\Desktop\tubo_troll.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\tubo_troll.stl');
stlPlot(vertices,faces,name)
%% filtering data
b=ones(1,10)/10; 
xfil = filter(b,1,x);
yfil = filter(b,1,y);
zfil = filter(b,1,z);
figure
plot3(x,y,z,'b',xfil(numel(b):end),yfil(numel(b):end),zfil(numel(b):end),'r');xlabel('x');ylabel('y')
                                                                            ;title('r filtered, b original')
x=xfil(numel(b):end);
y=yfil(numel(b):end);
z=zfil(numel(b):end);
%% get x y z
mystlPlot(vertices, faces, 0.3);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
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
%%
% con step=25 va
k=1;
step=10;
for i=10:step:numel(x)
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
plot3(x_plot(1,:),y_plot(1,:),z_plot(1,:),'Color',[0.5 0 i/numel(x)],'LineWidth',2);grid on;hold on 
k=k+1;
end

