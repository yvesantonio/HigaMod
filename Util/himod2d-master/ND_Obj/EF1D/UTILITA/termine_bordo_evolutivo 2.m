function r = termine_bordo_evolutivo(griglia,basi,dati_bordo)

  r = zeros(griglia.dim,1);
  if isnumeric(dati_bordo.r)
    valori = dati_bordo.r;
  else
    valori = feval(dati_bordo.r);
  end
  
  if dati_bordo.bc(1)~=1 
  % condizione di tipo Robin/Neumann: primo nodo
    r(1) = -valori(1);
  end
  if dati_bordo.bc(2)~=1  
  % condizione di tipo Robin/Neumann: secondo nodo
    r(end) = -valori(2);
  end
  
return
