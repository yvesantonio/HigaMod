[vertices,faces,normals,name] = stlRead('C:\Users\Leo\Desktop\my_model.stl');
hold on
mystlPlot(vertices,faces,0.3)
load('my_model_step5_from1_to_end_centerline2280punti.mat')
for i=250:1:numel(surf_point)
    coons_srf(i)=get_coons(surf_point{i},7,2);
   
    nrbplot(coons_srf(i),[10 10])
end
