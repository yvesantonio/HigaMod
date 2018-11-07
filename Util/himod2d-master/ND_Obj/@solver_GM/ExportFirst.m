function []=ExportFirst(obj)
            file=fopen('first','wt');
            fprintf(file,'%f\n%f\n%f\n%f\n',obj.first,obj.krpo,obj.full2d,obj.casetest);
            fclose(file);
        end