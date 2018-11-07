function MakeMoreImg(obj,v,varargin)
%
%  function MakeMoreImg(obj,v)
%  v deve essere un vettore che contiene gli istanti t_k che si vogliono
%  salvare come immagini.
%  Senza extra-argomenti salva in eps altrimenti salva nel formato dentro
%  varargin, per ora c'è solo png oltre ad eps.
%
if(~obj.solved)
    display('Necessario risolvere prima.')
    return
end
for i=1:length(v)
    obj.MakeImg(v(i),varargin);
end