function R = termine_bordo_costante(griglia,dati_bordo)
  
  R = sparse(griglia.dim,griglia.dim);
  
  if dati_bordo.bc(1)~=1 
  % condizione di tipo Robin/Neumann: primo nodo
      R(1,1) = -dati_bordo.gamma(1);
  end
  if dati_bordo.bc(2)~=1  
  % condizione di tipo Robin/Neumann: secondo nodo
      R(end,end) = -dati_bordo.gamma(2);
  end
  
return
