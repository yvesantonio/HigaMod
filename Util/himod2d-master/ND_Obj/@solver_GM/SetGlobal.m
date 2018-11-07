function []=SetGlobal(obj)
% function []=SetGlobal(obj)
%
% Questa funzione aggiorna le variabili global:
%       mu krpo larghezza Chi kappa Cin Cest w
% Utilizzare questa funzione ogniqualvolta si modifica la variabile
% dell'oggetto che rappresenta uno di questi valori
%            04/09/12
    global mu krpo larghezza Chi kappa Cin Cest
    mu=obj.mu;
    Chi=obj.Chi;
    kappa=obj.kappa;
    Cin=obj.Cin;
    Cest=obj.Cest;
    krpo=obj.krpo;
    larghezza=obj.larghezza;
end