%% stessi punti gradi differenti
clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0 ; 2 0; 2 2; 4 2; 5 1; 4 -1];% Pi=[x y]

figure; hold on; title('Changing the weight of a control points')
plot(P(:,1),P(:,2),'d')
p=2;
%creation of U considered P and p XDD
U=zeros(1,numel(P)/2+p+1); 
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
w=[0.8 0.9 1 1.1 1.2];
for i=1:5
 
N=nrbmak(P',U);
   N.coefs(4,4)=w(i);
nrbplot(N,1000)
end
grid on;xlabel('x'),ylabel('y')
%% surf
clear all
close all
clc
N=nrbtestsrf;
figure; hold on;grid on
nrbplot(N,[20 20],'colormap','jet');title('NURBS Surface')