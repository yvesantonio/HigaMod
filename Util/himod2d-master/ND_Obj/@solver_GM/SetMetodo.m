function []=SetMetodo(obj,str)
%   function []=SetMetodo(obj,str)
%
%   Setta il metodo di risoluzione come str che può assumere i valori:
%   krpo naive full2d
%   Se necessario ricalcola la mesh.
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem
            old=obj.full2d;
            switch str
                case 'full2d'
                    obj.full2d=1;
                    obj.krpo=0;
                case 'naive'
                    obj.full2d=0;
                    obj.krpo=0;
                case 'krpo'
                    obj.full2d=0;
                    obj.krpo=1;
                otherwise
                    error('there is no such method.');
            end
            obj.metodo=str;
            obj.SetGlobal;
            obj.ExportAll;
            obj.assembled=false;
end