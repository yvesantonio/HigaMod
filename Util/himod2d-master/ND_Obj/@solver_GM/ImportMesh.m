function []=ImportMesh(obj)
%    function []=ImportMesh(obj)
%    Salva le mesh coinvolte nella classe nelle variabili:
%    points2d,seg2d,tri2d
%    points1d,seg1d,tri1d
%    pointsfull2d,segfull2d,trifull2d
%    Bug:
%    Purtroppo salva nel Workspace principale delle variabili mesh.
%

if(~obj.full2d)
    obj.points2d=NaN;
    obj.seg2d=NaN;
    obj.tri2d=NaN;
    obj.points1d=NaN;
    obj.seg1d=NaN;
    obj.tri1d=NaN;
    
    [obj.points2d,obj.seg2d,obj.tri2d]=importfilemesh('mesh_2d.msh');
    [obj.points1d,obj.seg1d,obj.tri1d]=importfilemesh('mesh_1d.msh');
    
else
    obj.pointsfull2d=NaN;
    obj.segfull2d=NaN;
    obj.trifull2d=NaN;
    [obj.pointsfull2d,obj.segfull2d,obj.trifull2d]=importfilemesh('mesh_full2d.msh');
    
end
end