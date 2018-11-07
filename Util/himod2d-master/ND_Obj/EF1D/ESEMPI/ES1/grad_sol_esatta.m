function gradu = grad_sol_esatta(x)
  global beta

  C = (15+3/beta^2*(160*beta+31))/((beta+3)*exp(10*beta)-3);
  gradu = -beta*C*exp(beta*x) + 3/beta^2*(beta*x+1);
  
return
