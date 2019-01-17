clear all
close all
clc
s=@(x) sin(x);
c=@(x) cos(x);
t=0:0.01:2*pi;
x=s(t);
y=c(t);
z=s(4*t);
plot3(x,y,z,'b')