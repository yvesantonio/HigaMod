function []=ExportNeum(obj)
neum_per_2d=fopen('neum','wt');
fprintf(neum_per_2d,'%f\n%f\n',obj.val_neum_interfaccia,obj.val_neum_old);
fclose(neum_per_2d);
end