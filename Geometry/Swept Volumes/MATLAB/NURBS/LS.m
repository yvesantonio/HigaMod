f=@(x) (1-exp(-x));
dt=1/500; %frequenza campionamento scheda
t=0:dt:10; n=numel(t);
figure; plot(t,t,'b',t,f(t),'r');grid on;xlabel('t');ylabel('corrente');axis([0 10 0 1]);title('tangente in origine y=x')

I=f(t); %corrente ideale
d=randi([-1000 1000],1,n)*1.e-4; %distrubo che modella gli errori in genrale a media nulla
Ir=I+d; % corrnte reale = corrente e disturbo

figure; plot(t,f(t),'r',t,Ir,'b');grid on;xlabel('t');ylabel('corrente');axis([0 10 0 1]);title('tangente in origine y=x')
% leas sqaure del primo pezzo, dal grafico ideale hai che l'approssimazione
% è bella fino a 1/10 di tao
c=0:dt:0.05; %piu t finale è basso maggiore è l'approssimazione, ma piu sensibile al rumore :((

P=polyfit(c,Ir(1:numel(c)),1); %1 è lordine del polinomio che approssima
                                            %P restituisce i coefficienti
                                            %della retta
figure; plot(t,I,'r',t,(1-exp(-P(1)*t)));axis([0 10 0 1]);title('rosso reale, blue identificata') 
%ora piu prove 
ti=0; %istante di tempo iniziale per LS sempre 0
tf=0.05; %tempo finale
c=ti:dt:tf;
for i=1:1000
    Ir(i,:)=I+randi([-1000 1000],1,n)*1.e-5; %distrubo che modella gli errori in genrale a media nulla;
    Pi(i,:)=polyfit(c,Ir(i,1:numel(c)),1);
end
tao_stimato=sum(Pi(:,1))/numel(Pi(:,1));
figure; plot(t,I,'r',t,(1-exp(-tao_stimato*t)));axis([0 10 0 1]);title('rosso reale, blue identificata dopo 1000 prove') 