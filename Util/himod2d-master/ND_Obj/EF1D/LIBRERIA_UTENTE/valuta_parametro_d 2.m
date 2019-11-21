function d_ott = valuta_parametro_d(Pe)
% d_ott = valuta_parametro_d(Pe)
%
% Stima del parametro di viscosità numerica ottimale per
% stabilizzazione fortemente consistente.
%
% Si adotta una approssimazione della formula proposta in "Franca Frey
% Hughes: Stabilized finite element methods: application to the
% advective diffusive model, 1990"
%
% Pe: numero di Peclet di griglia
%
% N.B. a rigore la stima è riferita al caso rho=1

 if Pe<3
   d_ott = 0.5 * 1/3*Pe;
 else
   d_ott = 0.5;
 end

return

