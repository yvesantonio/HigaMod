clc
clear all
close all
%%
P(1,1,1,:)=[ 0 0 0 ]; P(2,1,1,:)=[ 1 0 0 ];  P(3,1,1,:)=[ 1 0 1 ];  P(4,1,1,:)=[ 2 0 1  ];  P(5,1,1,:)=[ 3 0 0 ];
P(1,1,2,:)=[ 0 0 1 ]; P(2,1,2,:)=[ 1 0 1 ];  P(3,1,2,:)=[ 1 0 2 ];  P(4,1,2,:)=[ 2 0 2  ];  P(5,1,2,:)=[ 3 0 1 ];
P(1,1,3,:)=[ 0 0 2 ]; P(2,1,3,:)=[ 1 0 2 ];  P(3,1,3,:)=[ 1 0 3 ];  P(4,1,3,:)=[ 2 0 3  ];  P(5,1,3,:)=[ 3 0 2 ];

P(1,2,1,:)=[ 0 1 -1 ];P(2,2,1,:)=[ 1 1 0 ];  P(3,2,1,:)=[ 1 1 0 ];  P(4,2,1,:)=[ 2 1 -1  ];  P(5,2,1,:)=[ 4 1 0 ];
P(1,2,2,:)=[ 0 1 1 ]; P(2,2,2,:)=[ 1 1 1 ];  P(3,2,2,:)=[ 2 1 1 ];  P(4,2,2,:)=[ 3 1 1  ];  P(5,2,2,:)=[ 4 1 1 ];
P(1,2,3,:)=[ 0 1 2 ]; P(2,2,3,:)=[ 1 1 2 ];  P(3,2,3,:)=[ 2 1 2 ];  P(4,2,3,:)=[ 3 1 2  ];  P(5,2,3,:)=[ 4 1 2 ];

figure
plot3(P(:,1,1,1), P(:,1,1,2),P(:,1,1,3),P(:,1,2,1), P(:,1,2,2),P(:,1,2,3),P(:,1,3,1), P(:,1,3,2),P(:,1,3,3));grid on;hold on
plot3(P(:,2,1,1), P(:,2,1,2),P(:,2,1,3),P(:,2,2,1), P(:,2,2,2),P(:,2,2,3),P(:,2,3,1), P(:,2,3,2),P(:,2,3,3));grid on
%%
U=0:0.1:1;
V=U;
W=U;
d=numel(U);


n=numel(P(:,1,1,1))-1;      % lunghezza di x
m=numel(P(1,:,1,1))-1;      % lunghezza di y
l=numel(P(1,1,:,1))-1;      % lunghezza di z
figure
%% direzione u

for nv=1:d
for nw=1:d
    C=zeros(3,d);
    for i=0:n

        bj=0;
        for j=0:m        
            bk=0;
            for k=0:l
                p(1:3,1)=P(i+1,j+1,k+1,:);
                bk=bk+bern(k,l,W(nw))*p;
            end
            bj=bj+bk*bern(j,m,V(nv));
        end
        C=C+[bern(i,n,U);bern(i,n,U);bern(i,n,U)].*bj;    
        
    end
    plot3(C(1,:),C(2,:),C(3,:),'r');grid on;xlabel('x');ylabel('y');hold on

end


end

%% direzione v

for nu=1:d
for nw=1:d
    C=zeros(3,d);
    for j=0:m

        bj=0;
        for i=0:n        
            bk=0;
            for k=0:l
                p(1:3,1)=P(i+1,j+1,k+1,:);
                bk=bk+bern(k,l,W(nw))*p;
            end
            bj=bj+bk*bern(i,n,U(nu));
        end
        C=C+[bern(j,m,U);bern(j,m,U);bern(j,m,U)].*bj;    
        
    end
    plot3(C(1,:),C(2,:),C(3,:),'g');grid on;xlabel('x');ylabel('y');hold on

end


end
%% direzione w

for nu=1:d
for nv=1:d
    C=zeros(3,d);
    for k=0:l

        bj=0;
        for i=0:n        
            bk=0;
            for j=0:m
                p(1:3,1)=P(i+1,j+1,k+1,:);
                bk=bk+bern(j,m,V(nv))*p;
            end
            bj=bj+bk*bern(i,n,U(nu));
        end
        C=C+[bern(k,l,W);bern(k,l,W);bern(k,l,W)].*bj;    
        
    end
    plot3(C(1,:),C(2,:),C(3,:),'b');grid on;xlabel('x');ylabel('y');hold on

end


end
%plot3(P(1:end,1:end,1:end,1),P(1:end,1:end,1:end,2),P(1:end,1:end,1:end,3),'b');grid on;xlabel('x');ylabel('y')

