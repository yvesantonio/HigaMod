function []=ExportGeom(obj)
            dati_geometrici=fopen('dati_geometrici','wt');
            fprintf(dati_geometrici,'%f\n%f\n%f\n%f\n%f\n%f\n%f\n',obj.lun2d,obj.lun1d,obj.larghezza,obj.spes,obj.dt,obj.timedipendent,obj.full2d);
            fclose(dati_geometrici);
        end