% Parte integrante della classe!!
obj.C = termine_trasporto(obj.griglia,obj.base,'coeff_a');

if(obj.krpo)
    obj.S = termine_reazione(obj.griglia,obj.base,'coeff_s');
end

dati_bordo = struct( ...
  'bc'    , [1 2] , ... % tipo delle condizioni 1->Dirichlet 2->Neumann
  'gamma' , [-1 -1] , ... % parametro gamma della condizione di Robin (-1 == non uso questo parametro)
  'r'     , [-1 0] ...  % dato di Neumann / Robin --> vedi guida (-1== non uso questo parametro)
);

%il primo grado di libertà è di Dirichlet
u_incognite = 2 : 1 : obj.griglia.dim;
u_note      = [1];


[R r] = termine_bordo(obj.griglia,obj.base,dati_bordo);

if(obj.timedipendent)
	A = obj.M/obj.dt+obj.theta*(obj.K+obj.C+R+obj.S);
	Z=-(1-obj.theta)*(obj.K+obj.C+R+obj.S)+obj.M/obj.dt;
	b = obj.bv+r+Z*obj.u1dold;
else	
	A = obj.K+obj.C+R+obj.S;
	b = obj.bv+r;
end
% Condizioni al contorno di Dirichlet: partizionamento del
% sistemali lineare

A11 = A(u_incognite,u_incognite);
A12 = A(u_incognite,u_note     );
A21 = A(u_note     ,u_incognite);
A22 = A(u_note     ,u_note     );
b1 = b(u_incognite);
x2 = [val_diric_interfaccia];  % valori noti di u al bordo

% risoluzione
x1 = A11\(b1-A12*x2);

% assemblo soluzione
u1dnew(u_incognite,1) = x1;
u1dnew(u_note     ,1) = x2;

% calcolo valore derivata all'interfaccia per problema 2D e grandezza
% d'arresto
aus=obj.val_neum_interfaccia;
obj.val_neum_interfaccia=(u1dnew(2)-u1dnew(1))/(obj.nodi(2)-obj.nodi(1));

errtest1dDominio=0;
for jei=1:length(obj.u1dkmeno1) %attenzione cambiato per lo stazionario possibile causa baco era obj.u1dold
	errtest1dDominio=errtest1dDominio+abs(obj.u1dkmeno1(jei)-u1dnew(jei))^2;
end

errtest1d=abs(aus-obj.val_neum_interfaccia);

%modifico la soluzione per comodità di plot
obj.u1d=zeros(length(u1dnew)*2,1);
for h=1:length(u1dnew)
    obj.u1d(2*h-1)=u1dnew(h);
    obj.u1d(2*h)=u1dnew(h);
end
obj.u1dkmeno1=u1dnew;
