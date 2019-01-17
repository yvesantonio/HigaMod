clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0 ; 1 -4; 2 0; 3 0; 4 2; 5 0; 6 0; 7 5; 8 0];% Pi=[x y]
p=3;
%creation of U considered P and p XDD
U=zeros(1,numel(P)/2+p+1); 
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
low=4; high=numel(U)-low;
u=0:0.001:0.999;
% i=6; % vettore parte da 1 non d a 0
% p=2; % basis degree
 p=3;
 figure;hold on
    for i=1:high
        if i==high
            span=numel(u);
        else
            span=numel(u)-1;
        end
    for k=1:span
        Nn(k)=Nip(i,u(k),p,U);
    end
    
    plot(u(1:span),Nn);hold on;title(['Basis ','for degree ',num2str(p)]); grid on;
    legend(['Basis for i=0'],['Basis for i=1'],['Basis for i=2'],['Basis for i=3'],['Basis for i=4'],['Basis for i=5'],['Basis for i=6'],...
        ['Basis for i=7'],['Basis for i=8'])
    end
figure; hold on; title('Changin one control point')
N=nrbmak(P',U);grid on;xlabel('x'),ylabel('y')
plot(P(:,1),P(:,2),'d')
nrbplot(N,1000)
P=[0 0 ; 1 -4; 2 0; 3 0; 2 4; 5 0; 6 0; 7 5; 8 0];
N=nrbmak(P',U);grid on;xlabel('x'),ylabel('y')
plot(P(:,1),P(:,2),'d')
nrbplot(N,1000)
%% stessi punti gradi differenti
clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0 ; 1 -3; 2 0; 3 0; 4 2; 5 0];% Pi=[x y]

figure; hold on; title('Changing the degree of the curve')
plot(P(:,1),P(:,2),'d')
for p=1:5
%creation of U considered P and p XDD
U=zeros(1,numel(P)/2+p+1); 
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
N=nrbmak(P',U);
nrbplot(N,1000)
end
grid on;xlabel('x'),ylabel('y');legend(['Control points'],['Curve of degree 1'],['Curve of degree 2'],['Curve of degree 3'],...
                                                                   ['Curve of degree 4'],['Curve of degree 5'])