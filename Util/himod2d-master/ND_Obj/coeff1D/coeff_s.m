%reazione
function s = coeff_s(x)

global mediaestremi Chi mu
global krpo
if(krpo)
    s=mediaestremi*Chi*mu*ones(size(x));
else
    s = zeros(size(x));
end
return
    
