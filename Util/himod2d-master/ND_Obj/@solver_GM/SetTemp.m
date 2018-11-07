function SetTemp(obj,t0,d,tf,theta,mtimedipendent)
%   function SetTemp(obj,t0,d,tf,theta,timedipendent)
%   
%   Setta l'istante iniziale,t0, il passo,d, il tempo finale,tf, e il theta
%   del theta-metodo, se mtimedipendent uguale a false lo memorizza e non fa
%	niente.
%   Aggiorna le variabili globali, mette assembled a falso ed aggiorna le
%   variabili condivise con freefem

obj.timedipendent=mtimedipendent;
if(mtimedipendent)
	obj.t_f=tf;
	obj.t_0=t0;
	obj.dt=d;
	obj.theta=theta;
end
obj.ExportAll;
obj.SetGlobal;
obj.assembled=false;
end
