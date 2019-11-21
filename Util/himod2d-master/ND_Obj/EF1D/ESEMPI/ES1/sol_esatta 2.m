function u = sol_esatta(x)
  global beta

  C = (15+3/beta^2*(160*beta+31))/((beta+3)*exp(10*beta)-3);
  u = C*(1-exp(beta*x)) + 3*x/(2*beta^2).*(beta*x+2);

return
