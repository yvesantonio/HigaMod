t=0:0.01:2;
x=2*sin(2*pi*t);
y=2*cos(2*pi*t);
z=2*t;
plot3(x,y,z,'b--','LineWidth',2)
plot_frenet(x,y,z,10)