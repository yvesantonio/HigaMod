clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% load data .dat , need to run the vmtk to get the .dat file
A=importdata('C:\Users\Leo\Desktop\cent_artcor_red.dat');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\artcor_red.stl');
mystlPlot(vertices,faces,1)
%% get x y z 
x=A.data(:,1);
y=A.data(:,2);
z=A.data(:,3);
r=A.data(:,4);
%%
mystlPlot(vertices, faces,1);
plot3(x,y,z,'b');grid on; xlabel('x'); ylabel('y')
% plot3(x(1),y(1),z(1),'rd');grid on; xlabel('x'); ylabel('y');
% plot3(x(end),y(end),z(end),'gd');grid on; xlabel('x'); ylabel('y');title('centerline b, inizio r, fine g')

%% prende quello che restituisce vmtk e ottiene us set di curve che vanno dall'inizo alle varie fini
jump(1)=1;
k=2;
for i=1:numel(x)-1
    if norm([x(i) y(i) z(i)]-[x(i+1) y(i+1) z(i+1)])>r(i)*2
        jump(k)=i;
        k=k+1;
    end
end
jump(end+1)=numel(x);

for i=1:numel(jump)-1
    C{i}=[x(jump(i+1):-1:jump(i)+1) y(jump(i+1):-1:jump(i)+1) z(jump(i+1):-1:jump(i)+1) r(jump(i+1):-1:jump(i)+1) ];
end
figure
% mystlPlot(vertices, faces, 0.1);
hold on
for i=1:numel(C)
    plot3(C{i}(:,1),C{i}(:,2),C{i}(:,3),'Color',[0 0 1],'LineWidth',2)
    plot3(C{i}(1,1),C{i}(1,2),C{i}(1,3),'ro')
    plot3(C{i}(end,1),C{i}(end,2),C{i}(end,3),'go')
end

%%

%%
% for i=[ 1 3 4]
%     [aaaaa N]=nrb_approx(C{i}(:,1),C{i}(:,2),C{i}(:,3),ceil(numel(C{i}(:,1))/20),3,1.e-3,10,10);
%     C{i}(:,1:3)=[nrbeval(N,linspace(0,1,numel(C{i}(:,1))))]';
% disp('ss')
% end

%%
St=getc(C); % crea la struttura ad albero segmentata

St=smoothcurve(St); %NURBS curve fitta i vari segmenti
mystlPlot_nofig(vertices, faces, 0.5);hold on
hold on
plotc(St)
%% plotta la tangente delle trene di frenet 
 figure
 p_fre(St.next2,2)
 axis('image');
 view([-135 35]);
 