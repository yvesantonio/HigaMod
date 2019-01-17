clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0; 1 1; 2 0; 3 1; 4 1; 5 -2; 6 0; 7 0; 8 1; 9 0]; % Pi=[x y]
% da P al nodo definito p a piacere
p=9;
U=[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1];
C=0;
u=0:0.001:1-0.001;
Cx=0;
Cy=0;
n=numel(P)/2;
for j=1:numel(u)
    Cx(j)=0;
    Cy(j)=0;
for i=1:n
    Cx(j)=Cx(j)+Nip(i,u(j),p,U)*P(i,1);    
    Cy(j)=Cy(j)+Nip(i,u(j),p,U)*P(i,2);    
end
end
figure; plot(P(:,1),P(:,2),'r',Cx,Cy,'b')
%%
clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0; 1 1; 2 0; 3 1; 4 1; 5 -2; 6 0; 7 0; 8 1; 9 0]; % Pi=[x y]
% da P al nodo definito p a piacere
p=2;
U=[0 0 0 1/8 2/8 3/8 4/8 5/8 6/8 7/8 1 1 1];
C=0;
u=0:0.001:1-0.001;
Cx=0;
Cy=0;
n=numel(P)/2;
for j=1:numel(u)
    Cx(j)=0;
    Cy(j)=0;
for i=1:n
    Cx(j)=Cx(j)+Nip(i,u(j),p,U)*P(i,1);    
    Cy(j)=Cy(j)+Nip(i,u(j),p,U)*P(i,2);    
end
end
figure; plot(P(:,1),P(:,2),'r',Cx,Cy,'b')
