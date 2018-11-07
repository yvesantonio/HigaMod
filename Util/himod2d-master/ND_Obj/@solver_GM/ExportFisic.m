function []=ExportFisic(obj)
            dati_fisici=fopen('dati_fisici','wt');
            fprintf(dati_fisici,'%f\n%f\n%f\n%f\n%f\n%f\n',obj.Cest,obj.Chi,obj.kappa,obj.Cin,obj.mu);
            fclose(dati_fisici);
        end