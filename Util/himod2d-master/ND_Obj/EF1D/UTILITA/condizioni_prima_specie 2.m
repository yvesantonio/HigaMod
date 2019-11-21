function x = condizioni_prima_specie(griglia,basi,dati_bordo)

  if isnumeric(dati_bordo.ub)
    valori = dati_bordo.ub;
  else
    valori = feval(dati_bordo.ub);
  end
  
  x = [];
  if dati_bordo.bc(1)==1  
  % condizione di tipo Dirichlet: primo nodo
    x = [x ; valori(1)];
  end
  if dati_bordo.bc(2)==1  
  % condizione di tipo Dirichlet: secondo nodo
    x = [x ; valori(2)];
  end
  
return
