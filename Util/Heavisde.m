% Heaviside Function

close all

Hd = @(x,d) 1./(1 + exp(-2 * x .* d));
func = @(x,y,d) (Hd(x+1,d) - Hd(x-1,d)).*(Hd(y+1,d) - Hd(y-1,d));

d = 100;
xx = linspace(-5,5,500);
yy = linspace(-5,5,500);
[X,Y] = meshgrid(xx,yy);

figure;
plot(xx,Hd(xx+1,d));

figure;
plot(xx,(Hd(xx+1,d) - Hd(xx-1,d)));

figure;
mesh(X,Y,func(X,Y,d));

% Circular Heaviside Function

Hc = @(x,y,a,d) 1./(1 + exp(-2 .* d .* (sqrt(x.^2 + y.^2) - a)));
func = @(x,y,a,d) Hc(x,y,-a,d) - Hc(x,y,a,d);

lag = 1;

figure;
surf(X,Y,func(X,Y,lag,d));