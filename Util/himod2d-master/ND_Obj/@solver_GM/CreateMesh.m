function []=CreateMesh(obj)
%   function []=CreateMesh(obj)
%   Chiama il file freefem ND_mesh.edp
%   si presuppone che tutto sia stato esportato
%   non modifica la classe
                !FreeFem++ -ne -nw ND_mesh0505.edp
end
