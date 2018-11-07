function [order,C]=compute_asimptotyc_order(uk,soglia)
%
% function [order,C]=compute_asimptotyc_order(uk,soglia)
% 
% dato un vettore uk e una soglia che indica di scartare gli elementi del vettore da 1 a soglia-1
% calcola order e C tali che:
%                          uk=CK^-order
%

k=soglia:length(uk);

[P]=polyfit(log10(k),log10(uk(k)'),1); % regressione lineare.

order=-P(1);
C=10^P(2);
