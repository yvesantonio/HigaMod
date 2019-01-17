clear all
close all
clc
% i span, nel libro i parte da 0, qua da 1
U=[0 0 0 0 3 5 7 9 9 9 9]/9; % more molteplicity more continuity :)
low=3; high=numel(U)-low;

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
    legend(['Basis for i=0'],['Basis for i=1'],['Basis for i=2'],['Basis for i=3'],['Basis for i=4'],['Basis for i=5'],['Basis for i=6'])
    end
    

%% derivatives
clear all
close all
clc
% i span, nel libro i parte da 0, qua da 1
U=[0 0  0 1  1 1  2 2 2  3 3 3  4 4 4  5 5 5]; % more molteplicity more continuity in the derivatives :)
low=3; high=numel(U)-low;

u=min(U):0.001:max(U);
i=7; % vettore parte da 1 non d a 0
p=3; % basis degree
figure
for k=1:numel(u)
    Nn(k)=Nip(i,u(k),p,U);
end

plot(u,Nn,'r');hold on;title(['Basis i=',num2str(i),' of degree ',num2str(p)])

for j=1:p
for k=1:numel(u)
        Nn(k)=DerNip(j,i,u(k),U,p);
end
figure
plot(u,Nn,'b');hold on;title(['Der ',num2str(j),'-th of the Basis i=',num2str(i),' of degree ',num2str(p)])
end





    
    
    
    
    
    