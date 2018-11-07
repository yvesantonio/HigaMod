function save_img(fig,name,varargin)
%
%  function save_img(fig,name,varargin)
%
%  senza extra-argomenti salva in png altrimenti salva nel formato dentro
%  name è parte del nome del file che viene salvato, il resto è data. 
% il tutto viene salvato dentro al cartella IMG
% fig è la figura da salvare
%  varargin, per ora c'è solo eps, ma funziona male oltre a png.
%
if(isempty(varargin))
    formato='png';
elseif(strcmp(varargin{1},'png'))
    formato='png';
elseif(strcmp(varargin{1},'eps'))
    formato='eps';
else
    display('in save_img: formato non riconosciuto')
    return
end
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];

titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    name ...
    );

if(strcmp(formato,'eps'))
    filename=strcat('./',titolo,'.eps');
    print(fig, '-depsc2', '-r300', filename);
else
    filename=strcat('./',titolo,'.png');
    print(fig, '-dpng', '-r300', filename);
end
end
