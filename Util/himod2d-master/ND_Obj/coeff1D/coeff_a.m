%trasporto
function a = coeff_a(x)

global krpo larghezza
if(krpo)
    global w;
    a=w*ones(size(x));
else
    a = 10/6*(larghezza^2)*ones(size(x));%//media
end
return
