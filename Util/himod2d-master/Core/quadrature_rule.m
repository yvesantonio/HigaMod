function [nod,wei]=quadrature_rule(a,b,nodes,weights)
%
% attualizza una formula di quadratura 1D da (0,1) all'intervallo (a,b)
% spostando i nodi e riscalando i pesi oppurtunamente
%

  h = b-a;
  nod = h*nodes+a;
  wei = h*weights;

return
