function [Map]=provideMap(haty)

invhaty= finverse(haty,'y');
Dnonhat=diff(haty,'x');
Jnonhat=diff(haty,'y');
D=subs(Dnonhat,'y',invhaty); % pay attention here y is becoming yhat
J=subs(Jnonhat,'y',invhaty); % here y is becoming yhat
Jac=diff(invhaty,'y'); %pay attention we derive in y but it should be yhav
a=sym('sqrt(1+invhaty_x^2)'); %%%!!!!!
invhaty_x=diff(invhaty,'x');
a=subs(a,'invhaty_x',invhaty_x);
a_up=subs(a,'y',1);
a_down=subs(a,'y',-1);

a_up=sym2fx(a_up);
a_down=sym2fx(a_down);
haty=sym2fxy(haty);
D=sym2fxy(D);
J=sym2fxy(J);
Jac=sym2fxy(Jac);
invhaty=sym2fxy(invhaty);
Map=struct('haty',haty,'invhaty',invhaty,'D',D,'J',J,'jac',Jac,'a_up',a_up,'a_down',a_down);

end