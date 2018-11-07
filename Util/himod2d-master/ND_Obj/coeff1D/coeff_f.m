%forzante
function f = coeff_f(x)
    global krpo mu Chi Cest larghezza
    if(krpo)
        f=2*mu*Chi*Cest*ones(size(x))/larghezza;
    else
        f = zeros(size(x));
    end
return
