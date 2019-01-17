clear all
close all
clc
%% point definition
P(1,1,1,:)=[ 0 0 0 ]; P(2,1,1,:)=[ 1 0 0 ];  P(3,1,1,:)=[ 1 0 1 ];  P(4,1,1,:)=[ 2 0 1  ];  P(5,1,1,:)=[ 3 0 0 ];
P(1,1,2,:)=[ 0 0 2 ]; P(2,1,2,:)=[ 1 0 2 ];  P(3,1,2,:)=[ 1 0 3 ];  P(4,1,2,:)=[ 2 0 3  ];  P(5,1,2,:)=[ 3 0 2 ];
P(1,1,3,:)=[ 0 0 3 ]; P(2,1,3,:)=[ 1 0 4 ];  P(3,1,3,:)=[ 1 0 5 ];  P(4,1,3,:)=[ 2 0 6  ];  P(5,1,3,:)=[ 3 0 5 ];

P(1,2,1,:)=[ 0 1 -1 ];P(2,2,1,:)=[ 1 1 0 ];  P(3,2,1,:)=[ 1 1 0 ];  P(4,2,1,:)=[ 2 1 -1  ];  P(5,2,1,:)=[ 4 1 0 ];
P(1,2,2,:)=[ 0 1 1 ]; P(2,2,2,:)=[ 1 1 1 ];  P(3,2,2,:)=[ 2 1 1 ];  P(4,2,2,:)=[ 3 1 1  ];  P(5,2,2,:)=[ 4 1 1 ];
P(1,2,3,:)=[ 0 1 2 ]; P(2,2,3,:)=[ 1 1 2 ];  P(3,2,3,:)=[ 2 1 2 ];  P(4,2,3,:)=[ 3 1 2  ];  P(5,2,3,:)=[ 4 1 2 ];

figure
plot3(P(:,1,1,1), P(:,1,1,2),P(:,1,1,3),P(:,1,2,1), P(:,1,2,2),P(:,1,2,3),P(:,1,3,1), P(:,1,3,2),P(:,1,3,3));grid on;hold on
plot3(P(:,2,1,1), P(:,2,1,2),P(:,2,1,3),P(:,2,2,1), P(:,2,2,2),P(:,2,2,3),P(:,2,3,1), P(:,2,3,2),P(:,2,3,3));grid on
%% automatic definition of equispaced knot 

pu=3; % degree of curve u
U=zeros(1,numel(P(:,1,1,1))+pu+1);
U(end-pu:end)=1;
mid=numel(U)-2*pu-2;
for i=1:mid
    U(i+pu+1)=1/(mid+1)*i;
end

pv=1;
V=zeros(1,numel(P(1,:,1,1))+pv+1);
V(end-pv:end)=1;
mid=numel(V)-2*pv-2;
for i=1:mid
    V(i+pv+1)=1/(mid+1)*i;
end

pw=2;
W=zeros(1,numel(P(1,1,:,1))+pw+1);
W(end-pw:end)=1;
mid=numel(W)-2*pw-2;
for i=1:mid
    W(i+pw+1)=1/(mid+1)*i;
end

u=[0:0.1:1-0.1 0.9999999999];
v=u;
w=u;
%%
for iw=1:numel(w)
    for iv=1:numel(v)
        for iu=1:numel(u)
            X(iu,iv,iw)=0;
            Y(iu,iv,iw)=0;
            Z(iu,iv,iw)=0; 
            for i=1:numel(P(:,1,1,1))
                n2=[0;0;0];
                    for j=1:numel(P(1,:,1,1))
                        n=0;
                        for k=1:numel(P(1,1,:,1))
                            p(1:3,1)=P(i,j,k,:);
                            n=n+Nip(k,w(iw),pw,W)*p;
                        end
                        n2=n2+n*Nip(j,v(iv),pv,V);
                    end
                X(iu,iv,iw)=X(iu,iv,iw)+n2(1)*Nip(i,u(iu),pu,U);
                Y(iu,iv,iw)=Y(iu,iv,iw)+n2(2)*Nip(i,u(iu),pu,U);
                Z(iu,iv,iw)=Z(iu,iv,iw)+n2(3)*Nip(i,u(iu),pu,U);
                
               
            end
        end
     end
end

%% u direction
figure
for iw=1:numel(w)
    for iv=1:numel(v)
        %plot3(  X(:,iv,iw),Y(:,iv,iw),Z(:,iv,iw),'r' );hold on;grid on
        plot3(  X(:,iv,iw),Y(:,iv,iw),Z(:,iv,iw),'Color',[1 iv/numel(v) iw/numel(w) ]  );hold on;grid on

    end
end
%% v direction
for iv=1:numel(v)
    for iu=1:numel(u)
        x(:)=X(iu,iv,:);
        y(:)=Y(iu,iv,:);
        z(:)=Z(iu,iv,:);
        %plot3(  x,y,z,'g'  );hold on;grid
        plot3(  x,y,z,'Color',[ iu/numel(u) iv/numel(v) 1 ]  );hold on;grid
    end
end
%% w direction
% for iu=1:numel(u)
%     for iw=1:numel(w)
%         x(:)=X(iu,:,iw);
%         y(:)=Y(iu,:,iw);
%         z(:)=Z(iu,:,iw);
%         plot3(  x,y,z,'b' ]  );hold on;grid
%         plot3(  x,y,z,'Color',[iu/numel(u) 1  iw/numel(w) ]  );hold on;grid
%     end
% end
%%
plot3(P(:,1,1,1), P(:,1,1,2),P(:,1,1,3),'k',P(:,1,2,1), P(:,1,2,2),P(:,1,2,3),'k',P(:,1,3,1), P(:,1,3,2),P(:,1,3,3),'k','LineWidth',3);grid on;hold on
plot3(P(:,2,1,1), P(:,2,1,2),P(:,2,1,3),'k',P(:,2,2,1), P(:,2,2,2),P(:,2,2,3),'k',P(:,2,3,1), P(:,2,3,2),P(:,2,3,3),'k','LineWidth',3);grid on