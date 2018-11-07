function MakeImg(obj,k,varargin)
%
%  function MakeImg(obj,k,varargin)
%
%  Salva l'istante t_k come immagine nella cartella IMG
%  senza extra-argomenti salva in ong altrimenti salva nel formato dentro
%  varargin, per ora c'è solo eps, ma funziona male oltre a png.
%
if(~obj.solved)
    display('Necessario risolvere prima.')
    return
end


if(isempty(varargin))
    formato='png';
elseif(strcmp(varargin{1},'png'))
    formato='png';
elseif(strcmp(varargin{1},'eps'))
    formato='eps';
else
    display('in MakeImg: formato non riconosciuto')
    return
end
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
if(obj.timedipendent)
    titolo=strcat(                                       ...
        date,'_',ora,                                            ...
        'caso',num2str(obj.casetest),                    ...
        obj.metodo,                                      ...
        '[',num2str(obj.t_0),',',num2str(obj.t_f),']',   ...
        'dt_',num2str(obj.dt),                           ...
        't_',num2str(k*obj.dt+obj.t_0)                   ...
        );
else
    titolo=strcat(                                       ...
        date,'_',ora,                                            ...
        'caso',num2str(obj.casetest),                    ...
        obj.metodo,                                      ...
        'steady'                                           ...
        );
end
obj.currentfig=figure;
obj.PlotSingleTimeStep(k);

if(strcmp(formato,'eps'))
    filename=strcat(pwd,'/IMG/',titolo,'.eps');
    print(obj.currentfig, '-depsc2', '-r300', filename);
else
    filename=strcat(pwd,'/IMG/',titolo,'.png');
    print(obj.currentfig, '-dpng', '-r300', filename);
end
close(obj.currentfig);
end