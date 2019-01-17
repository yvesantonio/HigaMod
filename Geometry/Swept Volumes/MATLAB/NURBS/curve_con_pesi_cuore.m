clear all
close all
clc
%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0; 2 2; 4 0; 1 1 ];% Pi=[x y]
p=2;
%creation of U considered P and p XDD
U=zeros(1,numel(P)/2+p+1); % U=[0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1 1 1];
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
C=0;
u=[0:0.01:1-0.01 0.999];
n=numel(P)/2;
for j=1:numel(u)
    Cx(j)=0;
    Cy(j)=0;
        for i=1:n
            Cx(j)=Cx(j)+Nip(i,u(j),p,U)*P(i,1);    
            Cy(j)=Cy(j)+Nip(i,u(j),p,U)*P(i,2);    
        end
end


figure
plot(P(:,1),P(:,2),'bo',Cx,Cy,'r','MarkerSize',10);

%%
clear all
close all
clc
figure 
hold on

%%U=[0 0  0 0.25 0.25  0.5 0.75 1 1 1];
P=[0 0 ; 2 0; 2 2; ; 4 2; 5 1; 4 -1];% Pi=[x y]
W=[ 1    1           1      1     1    1 ];% weight
p=3;
wee=[0 0.5 1 2 10];
%creation of U considered P and p XDD
for iw=1:numel(wee)
    W(4)=wee(iw);
U=zeros(1,numel(P)/2+p+1); 
U(end-p:end)=1;
mid=numel(U)-2*p-2;
for i=1:mid
    U(i+p+1)=1/(mid+1)*i;
end
C=0;
u=[0:0.005:1-0.005 0.99999];
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
plot(P(:,1),P(:,2),'rd',Cx/w,Cy/w,Cx(130)/w,Cy(130)/w,'d');
end
