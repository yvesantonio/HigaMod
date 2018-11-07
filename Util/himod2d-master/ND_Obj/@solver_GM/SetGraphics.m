function SetGraphics(obj,treD,frontale,az,el,varargin)
%   function SetGraphics(obj,treD,frontale,az,el,cinf,csup)
%
%   treD=1 se si vuole un grafico 3D
%   frontale vale solo se treD=1 per guardare frontalmente la soluzione
%   az,el rispettivamente angolo orizzontale e verticale
%   cinf, csup : range scala di colore.
%
%   function SetGraphics(obj,treD,frontale,az,el,stringa)
%
%   se in stringa si mette auto il range della scala di colore viene
%   settato in base alla soluzione corrente
%
%   function SetGraphics(obj,treD,frontale,az,el)
%
%   vengono mantenuti i valori correnti di cinf csup
obj.treD=treD;
obj.frontale=frontale;
obj.az=az;   %Angolo visuale orizzontale
obj.el=el;     %Angolo visuale verticale
if(size(varargin,2)==1)
    if(strcmp(varargin{1},'auto'))
        if(obj.solved)
             obj.cinf=+inf;
             obj.csup=-inf;
            for i=1:size(obj.Solution,1)
                if(obj.full2d)
                    obj.cinf=min(min(obj.Solution{i,1}),obj.cinf);
                    obj.csup=max(max(obj.Solution{i,1}),obj.csup);
                else
                    obj.cinf=min(min(min(obj.Solution{i,1}),obj.cinf),min(obj.Solution{i,2}));
                    obj.csup=max(max(max(obj.Solution{i,1}),obj.csup),max(obj.Solution{i,2}));
                end
            end
        else
            display('cinf,csup invariati, per usare la funzione auto prima risolvere il problema');
        end
    else
        display('errore fornire due valori uno per cinf e uno csup')
    end
elseif(size(varargin,2)==2)
    if(isnumeric(varargin{1}) && isnumeric(varargin{2}))
        obj.cinf=varargin{1};
        obj.csup=varargin{2};
       
    else
        display('gli argomenti aggiuntivi possono essere o auto o due valori numerici')
    end
elseif(size(varargin,2)>2)
    display('gli argomenti aggiuntivi possono essere o auto o due valori numerici')
end

obj.levels=linspace(obj.cinf,obj.csup,20);

end