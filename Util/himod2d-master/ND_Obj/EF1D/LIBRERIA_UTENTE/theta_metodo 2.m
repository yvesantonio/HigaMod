function [tempi,u,varargout] = theta_metodo(griglia,base,M,Ai, ...
                                            dati_bordo,fvol,u0, ...
                                            theta,Tiniz,Dt,Tfin)
% [tempi,u] = theta_metodo(griglia,base,M,Ai, ...
%                          dati_bordo,fvol,u0, ...
%                          theta,Tiniz,Dt,Tfin)
%
% Evoluzione temporale per problemi parabolici mediante theta-metodo
%
% Possono dipendere dal tempo i valori delle condizioni al bordo e il 
% valore del termine noto. Fa eccezione il parametro 'gamma' per 
% condizioni di tipo Robin, che deve essere un valore costante. 
% Non possono invece dipendere dal tempo il tipo delle condizioni al 
% bordo ed i coefficienti.
%
% griglia: struttura che descrive la griglia di calcolo, ottenuta
%          con 'struttura_griglia'
% base: struttura con le informazioni circa le funzioni di forma,
%       ottenuta con 'struttura_base'
% M: matrice di massa per il problema (assunta costante)
% Ai: matrice associata alla semidiscretizzazione sapaziale del problema
%     In tale matrice devono essere presenti tutti i termini associati
%     agli operatori di derivazione spaziale ad eccezione dell'eventuale
%     termine dato da condizioni al bordo di III specie (il termine R
%     restituito dalla funzione 'termine_bordo')
% dati_bordo: struttura con le informazioni sulle condizioni al bordo.
%             N.B. tale struttura differisce dalla analoga struttura 
%             utilizzata per problemi stazionari, in quanto deve essere
%             possibile specificare condizioni dipendenti dal tempo.
%             Essa si compone dei tre campi seguenti:
%             bc: vettore di due elementi corrispondenti agli estremi 
%                 dell'intervallo di definizione del problema. Ciascun 
%                 elementi vale 1 se la corrispondente condizione e' 
%                 di Dirichlet, 2 se di Neuman, 3 per il caso Robin.
%             gamma: vettore di due elementi con i coefficienti gamma 
%                    della condizione di Robin ai due estremi. N.B. e'
%                    necessario assegnare un valore anche in presenza
%                    condizioni di Diriclhet. Tale valore verra' poi 
%                    ignorato.
%             r: stringa contenete il nome della funzione che valuta il
%                dato al bordo. Tale funzione non ha parametri di
%                ingresso e restituisce un vettore di due elementi con 
%                il valore del dato r per condizioni al bordo di Neumann
%                e Robin. N.B. e' necessario assegnare un valore anche 
%                in presenza condizioni di Diriclhet o di Neumann. Tale
%                valore verra' poi ignorato.
%                Se il dato r e' costante, e' possibile specificarne 
%                direttamente il valore come vettore di due elementi.
%             ub: stringa contenete il nome della funzione che valuta il
%                 dato al bordo. Tale funzione non ha parametri di
%                 ingresso e restituisce un vettore di due elementi con 
%                 il valore del dato r per condizioni al bordo di 
%                 Dirichlet. N.B. e' necessario assegnare un valore 
%                 anche in presenza condizioni di Neumann o di Robin. 
%                 Tale valore verra' poi ignorato.
%                 Se il dato ub e' costante, e' possibile specificarne 
%                 direttamente il valore come vettore di due elementi.
% fvol: stringa con il nome della funzione che valuta il termine noto.
%       Tale funzione riceve in ingresso un vettore di coordinate x e
%       restituisce il valore del dato f come vettore della stessa
%       dimensione di x.
%       N.B. e' anche possibile utilizzare un puntatore a funzione
% u0: dato iniziale (come vettore colonna)
% theta: parametro theta in [0 1]
% Tiniz: istante iniziale
% Dt: passo temporale
% Tfin: istante finale
% 
% tempi: istanti temporali per i quali viene calcolata X (vettore riga)
% u: soluzione (ciascuna colonna corrisponde ad un istante temporale)
%
% [tempi,X,B] = theta_metodo(griglia,base, ...
%                            M,Ai,dati_bordo,dati_volume,x0, ...
%                            theta,Tiniz,Dt,Tfin)
% B: matrice di iterazione. Condizione necessaria e sufficiente perche'
%    il metodo sia stabile e' che il raggio spettrale di B sia minore di
%    1. Condizione sufficiente per la stabilita' e' norm(B)<1.
%
% N.B. per specificare dati dipendenti dal tempo e' disponibile la 
% variabile globale "TEMPO" definita nella funzione.
 
 % Controllo sui dati
 if(~isfield(dati_bordo,'ub'))
    error(['Per problemi evolutivi e'' necessario aggiungere ', ...
           'un campo ''ub'' alla struttura ''dati_bordo''']);
 end
 if or(~isnumeric(dati_bordo.gamma(1)), ...
       ~isnumeric(dati_bordo.gamma(2)))
   error('Attenzione: il parametro ''gamma'' deve essere una costante');
 end
 
 global TEMPO
 
 % Ispecie indica la presenza di almeno una condizione di I specie
 Ispecie = 0;
 
 if and(dati_bordo.bc(1)==1 , dati_bordo.bc(2)==1)
    u_note = [1 griglia.dim];
    u_incognite = [2 : 1 : griglia.dim-1];
    Ispecie = 1; 
 elseif and(dati_bordo.bc(1)==1 , dati_bordo.bc(2)~=1)
    u_note = [1];
    u_incognite = [2 : 1 : griglia.dim];
    Ispecie = 1;
 elseif and(dati_bordo.bc(1)~=1 , dati_bordo.bc(2)==1)
    u_note = [griglia.dim];
    u_incognite = [1 : 1 : griglia.dim-1];
    Ispecie = 1;
 else
    u_note = [];
    u_incognite = [1 : 1 : griglia.dim];
 end
 
 % Il termine di bordo e' spezzato in due contributi: il seguente,
 % costante, che viene inserito nella matrice del sistema, ed
 % un termine noto che puo' essere variabile.
 R = termine_bordo_costante(griglia,dati_bordo);
 A = Ai + R;
 
 B = M + Dt*theta*A;
 Z = M - Dt*(1-theta)*A;
 
 B11 = B(u_incognite,u_incognite);
 B12 = B(u_incognite,u_note);
 
 Z1 = Z(u_incognite,:);
 
 % essendo la matrice costante e' conveniente fattorizzarla
 % una volta per tutte
 [L U P] = lu(B11);
 
 % condizioni iniziali
 TEMPO = Tiniz;
 bvol = termine_noto(griglia,base,fvol);
 r = termine_bordo_evolutivo(griglia,base,dati_bordo);
    
 b = bvol+r;
 b1v = b(u_incognite);
 uk = u0;
 
 tempi(1) = TEMPO;
 u(:,1) = uk;
 
 iter = 0;
 for TEMPO=Tiniz+Dt:Dt:Tfin
    
    iter = iter+1;  % contatore
    
    % valutazione dei coefficienti
    bvol = termine_noto(griglia,base,fvol);
    r = termine_bordo_evolutivo(griglia,base,dati_bordo);
    
    b = bvol+r;
    b1 = b(u_incognite);
    
    if Ispecie
      u2 = condizioni_prima_specie(griglia,base,dati_bordo);
    
      % risoluzione
      y1 = L\(P*(Z1*uk + theta*b1 + (1-theta)*b1v - B12*u2));
      u1 = U\y1;
      uk(u_incognite,1) = u1;
      uk(u_note,1) = u2;
      
    else
    
      % risoluzione
      y1 = L\(P*(Z1*uk + theta*b1 + (1-theta)*b1v));
      uk = U\y1;
      
    end
    
    tempi(iter+1) = TEMPO;
    u(:,iter+1) = uk;
    
    b1v = b1;
    
 end
 
 if nargout==3
     % richiesta la matrice di iterazione
     varargout = {inv(B11)*Z(u_incognite,u_incognite)};
 end
 
return
    
