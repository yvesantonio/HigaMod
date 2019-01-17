clc
clear all
close all
U=0:0.01:1;
V=U;
d=numel(U);
X=zeros(d,d);
Y=zeros(d,d);
Z=zeros(d,d);

P(:,: ,1)=[0 0 0;1 0 0; 2 0 0;3 0 0];
P(:,: ,2)=[0 1 0;1 1 1; 2 1 -1;3 1 0];
P(:,: ,3)=[0 2 0;1 2 1; 2 2 2;3 2 3];
n=numel(P(1,:,1))-1;
m=numel(P(1,1,:))-1;
figure
plot3(P(:,1,1),P(:,2,1),P(:,3,1),P(:,1,2),P(:,2,2),P(:,3,2),P(:,1,3),P(:,2,3),P(:,3,3));grid on;xlabel('x');ylabel('y')
figure
for k=1:d
    u0=V(k);
    C=zeros(3,d);
for i=0:m 
    b0=0;
    for j=0:n
        b0=b0+bern(j,m,u0)*P(j+1,:,i+1)';
    end
    o=bern(i,n,U);
    C=C+b0.*[o;o;o];
end
plot3(C(1,:),C(2,:),C(3,:));grid on;xlabel('x');ylabel('y');hold on
end
