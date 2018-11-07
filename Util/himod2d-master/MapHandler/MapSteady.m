function maps=MapSteady(map,t)

haty=@(x,y) map.haty(t,x,y);
invhaty=@(x,y) map.invhaty(t,x,y);
D=@(x,y) map.D(t,x,y); % pay attention here y is becoming yhat
J=@(x,y) map.J(t,x,y); % here y is becoming yhat
Jac=@(x,y) map.Jac(t,x,y); %pay attention we derive in y but it should be yhav
a_up=@(x,y) map.a_up(t,x,y);
a_down=@(x,y) map.a_down(t,x,y);

maps=struct('haty',haty,'invhaty',invhaty,'D',D,'J',J,'jac',Jac,'a_up',a_up,'a_down',a_down);

end