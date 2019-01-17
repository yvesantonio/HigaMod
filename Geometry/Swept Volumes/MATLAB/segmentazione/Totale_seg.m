clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% load data .dat , need to run the vmtk to get the .dat file
A=importdata('C:\Users\Leo\Desktop\cent_model_complete.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\model.stl');
stlPlot(vertices,faces,name)
%% get x y z 
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
r=A.data(:,4);
%%
mystlPlot(vertices, faces,0.1);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

%%
jump(1)=1;
k=2;
for i=1:numel(x)-1
    if norm([x(i) y(i) z(i)]-[x(i+1) y(i+1) z(i+1)])>r(i)
        jump(k)=i;
        k=k+1;
    end
end
jump(end+1)=numel(x);

for i=1:numel(jump)-1
    C{i}=[x(jump(i+1):-1:jump(i)+1) y(jump(i+1):-1:jump(i)+1) z(jump(i+1):-1:jump(i)+1) r(jump(i+1):-1:jump(i)+1) ];
end

mystlPlot(vertices, faces, 0.1);hold on
for i=1:numel(C)
    plot3(C{i}(:,1),C{i}(:,2),C{i}(:,3),'Color',[1-i/numel(C) 0 i/numel(C)])
    plot3(C{i}(1,1),C{i}(1,2),C{i}(1,3),'rd')
    plot3(C{i}(end,1),C{i}(end,2),C{i}(end,3),'bd')
end
pause(1)
%%
St=getc(C);
%%
mystlPlot(vertices, faces, 0.5);hold on
hold on
plotc(St)
pause(1)
%%
f=faces;
v=vertices;
[f v]=closefaces(f,v);
mystlPlot(v,f,1)
disp('pdd')
[St f v]=create_section(St,v,f);
f=faces;
v=vertices;


[St f v]=add_mesh(St,f,v);

mystlPlot(v, f, 0.3);hold on
plotlimit(St,v)
plotc(St)
%plotsection(St)
pause(1)
%% segmentation
St.tuboexist=1;
St.start=1;
St=segmentazione(St,f,v);

%% plot segment
figure;hold on
plot_segment(St,v)