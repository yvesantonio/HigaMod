function MakeVideo(obj)
%    function MakeVideo(obj)
%    
%     Fa una video con la soluzione corrente e il grafico corrente.
%     viene salvato in pwd/VIDEO/ con un nome che è una stringa
%     automatica così composta:
%     date,'caso',obj.casetest,obj.metodo,[,obj.t_0,obj.t_f,]obj.dt                            ...
%
if(~obj.solved)
    display('Necessario risolvere prima.')
    return
end
orario=clock();
ora=[num2str(orario(4)),':',num2str(orario(5))];
titolo=strcat(                                       ...
    date,'_',ora,                                            ...
    'caso',num2str(obj.casetest),                    ...
    obj.metodo,                                      ...
    '[',num2str(obj.t_0),',',num2str(obj.t_f),']',   ...
    'dt_',num2str(obj.dt)                            ...
    );
n=size(obj.Solution,1);
obj.currentfig=figure;
nomefilmato=strcat(pwd,'/VIDEO/',titolo,'.avi');
filmatoavi = avifile(nomefilmato,'fps',5); %fps basso=filmato lento.
for k=1:n
    obj.PlotSingleTimeStep(k)
    F= getframe(obj.currentfig);
    filmatoavi=addframe(filmatoavi,F);
    close(obj.currentfig)
end
filmatoavi=close(filmatoavi);
end