PP=[0 0;
   1 2;
   3 1;
   2 0;
   0 0];
P1=PP(:,2);
u=0:0.01:1;
% calcolo B
n=numel(P1)-1;
figure; hold on
    for i=0:n
        B(i+1,:)=factorial(n)/factorial(i)/factorial(n-i)*u.^i.*(1-u).^(n-i);

    end
figure
plot(0:0.01:1,(0:0.01:1)*2,0,0,'rd',1,2,'rd',3,1,'rd',2,0,'rd',PP(1,1)*B(1,:)+PP(2,1)*B(2,:)+PP(3,1)*B(3,:)+PP(4,1)*B(4,:),PP(1,2)*B(1,:)+PP(2,2)*B(2,:)+PP(3,2)*B(3,:)+PP(4,2)*B(4,:))

%%
% loop?
clc
clear all
close all
P1=[0 0; 1 1  ; 2 -3 ;  5 4 ; 6 0 ; 4 -1 ];

u=0:0.01:1;
% calcolo B
n=numel(P1(:,1))-1;
figure; hold on; grid on;xlabel('u');ylabel('Basis value');title('Basis function for 5-th degree Bèzier Curve')
    for i=0:n
        B1(i+1,:)=factorial(n)/factorial(i)/factorial(n-i)*u.^i.*(1-u).^(n-i);
                plot(u,B1(i+1,:));legend(['Basis i=0'],['Basis i=1'],['Basis i=2'],...
                                        ['Basis i=3'],['Basis i=4'],['Basis i=5'])
    end
    c1=zeros(3,numel(u));
% Y=f(X)
for r=1:2 % x y z
for k=1:n+1 % tutti i punti
    c1(r,:)=c1(r,:)+P1(k,r)*B1(k,:);
end
end
figure

%plot(P(:,1),P(:,2),C(1,:),C(2,:))
plot(c1(1,:),c1(2,:));grid on;xlabel('x');ylabel('y');zlabel('z');hold on
plot(c1(1,70),c1(2,70),'rd')
plot(P1(:,1),P1(:,2),'rd')
P1=[0 0; 1 1  ; 2 -3 ;  4 5 ; 6 0 ; 4 -1 ];
    c1=zeros(3,numel(u));
for r=1:2 % x y z
for k=1:n+1 % tutti i punti
    c1(r,:)=c1(r,:)+P1(k,r)*B1(k,:);
end
end
plot(c1(1,:),c1(2,:));grid on;xlabel('x');ylabel('y');zlabel('z');hold on
plot(c1(1,70),c1(2,70),'rd')
plot(P1(:,1),P1(:,2),'rd')
% plot3(c2(1,:),c2(2,:),c2(3,:));hold on
% plot3(c3(1,:),c3(2,:),c3(3,:));hold on
% generazione p int
% pi1=c1(:,50);
% pi2=c2(:,50);
% pi3=c3(:,50);
% PY=polyfit([pi1(1) pi2(1) pi3(1)],[pi1(2) pi2(2) pi3(2)],2);
% PZ=polyfit([pi1(1) pi2(1) pi3(1)],[pi1(3) pi2(3) pi3(3)],2);
