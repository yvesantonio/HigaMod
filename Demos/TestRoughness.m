w = 0.2;
phi = 2*pi;
theta1 = 0;
M1 = 20;
k1 = 1;

g = @(s) cos((w * pi * s).^3) .* exp(-( (phi * s - theta1)/(M1 * k1)).^2);
f = @(x,y) (g(x) .* g(y));

xx = linspace(-2*pi,2*pi,500);
yy = linspace(-2*pi,2*pi,500);
[X,Y] = meshgrid(xx,yy);
ff = f(X,Y);

figure;
mesh(X,Y,ff);