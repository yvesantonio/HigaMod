kappa =2;
larghezza=1;
delta=0.5;
f = @(y) 20*(kappa+2)/(kappa+2*(1-delta))*(1-delta*((y-0.5)*2/larghezza).^kappa); 
quad(f,0,1)