clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\ArtCor.stl');
stlPlot(vertices,faces,name); xlabel('x');ylabel('y');zlabel('z')
k=1;
for i=1:numel(faces)/3
    if (vertices(faces(i,:),1)<-4.5 & vertices(faces(i,:),1)>-30) & ...
       (vertices(faces(i,:),2)<20 & vertices(faces(i,:),2)>-40) & ...
       (vertices(faces(i,:),3)<20 & vertices(faces(i,:),3)>-45) 
   new_f(k,:)=[faces(i,:)];
   k=k+1;
    end
end

figure
stlPlot(vertices,new_f,name)

stlWrite('C:\Users\Leo\Desktop\ArtCor_red.stl',new_f,vertices)
