function base = aggiunta_lapPHI(base)
% base = aggiunta_lapPHI(base)
%
% Funzione per l'aggiunta alla struttura base del campo 'lapPHI',
% contenente la derivata seconda (laplaciano) delle funzioni di forma.
%
% NOTA: tipicamente le formulazioni agli elementi finiti non
% necessitane le derivate seconde delle funzioni di forma, tali
% derivate sono però utilizzate in casi particolari, quali ad esempio
% le formulazioni stabilizzate fortemente consistenti.
%
% base: struttura ottenuta con la funzione 'struttura_base' alla quale
%       si desideri aggiungere il campo 'lapPHI'

 [dummy, dummy, base.lapPHI] = LagrPoli(base.csinodi, base.csiGauss);

return

