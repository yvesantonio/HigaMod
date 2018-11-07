function []=SaveTimeStep(obj,num)
%    function []=SaveTimeStep(obj,num)
%
%   Questa funzione salva nella variabile solution la soluzione al passo i
%   i=1 si intende t=0;
%

if(obj.timedipendent)
    if(~obj.full2d)
        obj.Solution{num,1}=importfiledata('u2dold');
        obj.Solution{num,2}=obj.u1d;
        if(obj.krpo)
            obj.Import_profile(num);
        end
    else
        obj.Solution{num,1}=importfiledata('u2dold');
        obj.Solution{num,2}=NaN;
    end
else
    if(~obj.full2d)
        obj.Solution{num,1}=importfiledata('u2dkmeno1');
        obj.Solution{num,2}=obj.u1d;
        if(obj.krpo)
            obj.Import_profile(num);
        end
    else
        obj.Solution{num,1}=importfiledata('u2dkmeno1');
        obj.Solution{num,2}=NaN;
    end
end
end