% clear all
% close all
% clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0 ; 2 0; 2 2; 4 2; 5 1; 3 -1]';
x=P(1,:);
y=P(2,:);
z=y*0;% Pi=[x y]
p=2;
%creation of U considered P and p XDD
U=zeros(1,numel(x)+p+1); % U=[0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1 1 1];
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
C=0;
u=linspace(0,0.999,10000);
n=numel(P)/3;
for j=1:numel(u)
    Cx(j)=0;
    Cy(j)=0;
    Cz(j)=0;
        for i=1:n
            Cx(j)=Cx(j)+Nip(i,u(j),p,U)*P(i,1);    
            Cy(j)=Cy(j)+Nip(i,u(j),p,U)*P(i,2);  
            Cz(j)=Cz(j)+Nip(i,u(j),p,U)*P(i,3);
        end
end


figure
plot3(P(:,1),P(:,2),P(:,3),'b',Cx,Cy,Cz,'r','MarkerSize',10);grid on

%%
clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0; 2 2; 4 0; 5 10; 1 1];% Pi=[x y]
W=numel(P)/2;% weight
p=1;
%creation of U considered P and p XDD
U=zeros(1,numel(P)/2+p+1); 
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
C=0;
u=[0:0.001:1-0.001 0.9999];
n=numel(P)/2;
for j=1:numel(u)
    Cx(j)=0;
    Cy(j)=0;
        for i=1:n
            Cx(j)=Cx(j)+Nip(i,u(j),p,U)*P(i,1)*W(i);    
            Cy(j)=Cy(j)+Nip(i,u(j),p,U)*P(i,2)*W(i);    
        end
    w=0;
    for i=1:n
        w=w+Nip(i,u(j),p,U)*W(i);
    end
    Cx(j)=Cx(j)/w;
    Cy(j)=Cy(j)/w;
end

figure
plot(P(:,1),P(:,2),'r',Cx/w,Cy/w,'b');
