clc
clear all
close all
P(:,1,1,:)=[0 0 0; 1 0 1; 2 0 2; 3 0 0];    P(:,2,1,:)=[0 1 0; 1 1 1; 2 1 2; 3 1 0];
P(:,1,2,:)=[0 0 1; 1 0 2; 2 0 3; 3 0 1];    P(:,2,2,:)=[0 1 1; 1 1 2; 2 1 3; 3 1 1];
%%
U=0:0.1:1;
V=U;
W=U;
d=numel(U);


n=numel(P(:,1,1,1))-1;      % lunghezza di x
m=numel(P(1,:,1,1))-1;      % lunghezza di y
l=numel(P(1,1,:,1))-1;      % lunghezza di z
figure
plot3(P(:,1,1,1),P(:,1,1,2),P(:,1,1,3),'g',P(:,2,1,1),P(:,2,1,2),P(:,2,1,3),'g',P(:,1,2,1),P(:,1,2,2),P(:,1,2,3),'g',P(:,2,2,1),P(:,2,2,2),P(:,2,2,3),'g');grid on;xlabel('x');ylabel('y')
figure
%%
for iw=1:numel(w)
    for iv=1:numel(v)
        for iu=1:numel(u)
            X(iu,iv,iw)=0;
            Y(iu,iv,iw)=0;
            Z(iu,iv,iw)=0;
            for i=1:numel(P(:,1,1,1))
                for j=1:numel(P(1,:,1,1))
                for k=1:numel(P(1,1,:,1))
                        X(iu,iv,iw)=X(iu,iv,iw)+Nip(j,u(iu),pu,U)*Nip(i,v(iv),pv,V)*Nip(k,w(iw),pw,W)*P(i,j,k,1);
                        Y(iu,iv,iw)=Y(iu,iv,iw)+Nip(j,u(iu),pu,U)*Nip(i,v(iv),pv,V)*Nip(k,w(iw),pw,W)*P(i,j,k,2);
                        Z(iu,iv,iw)=Z(iu,iv,iw)+Nip(j,u(iu),pu,U)*Nip(i,v(iv),pv,V)*Nip(k,w(iw),pw,W)*P(i,j,k,3);
                end
            end
            end
        end
    end
end
