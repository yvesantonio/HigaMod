function SetGeom(obj,mlun2d,mlun1d,mlarghezza,mspes)
%dati geometrici
        obj.lun2d=mlun2d;                    % attenzione deve coincidere nel file ND_mesh.edp
        obj.lun1d=mlun1d;                    % attenzione deve coincidere nel file ND_mesh.edp
        obj.larghezza=mlarghezza;            %in millimetri
        obj.spes=mspes;
        obj.ExportAll;
        %obj.CreateMesh;
        %obj.ImportMesh;
        obj.assembled=false;
end
