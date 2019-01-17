hold on
mystlPlot(vertices,faces,0.3)
[vertices,faces,normals,name] = stlRead('sphere_ascii.stl');
stlPlot(vertices,faces,name);

load('section_model.mat')
for i=1:numel(surf_point)
    coons_srf(i)=get_coons(surf_point{i},6,2);
   
    nrbplot(coons_srf(i),[20 20])
end
