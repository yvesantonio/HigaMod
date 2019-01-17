clear all
close all
clc
%% definizione superficie tramite punti versp il basso X; verso destra Y; P(x,y,z)
P(1,1,:)=[0 0 1];   P(1,2,:)=[0 1 1];  P(1,3,:)=[0 2 0.8];   P(1,4,:)=[0 2 0.8]; P(1,5,:)=[0 2 0.8];          P(1,6,:)=[0 3 1];     P(1,7,:)=[0 4 1];
P(2,1,:)=[1 0 0.8];   P(2,2,:)=[1 1 0];   P(2,3,:)=[1 2 0.7]; P(2,4,:)=[1 2 0.7];P(2,5,:)=[1 2 0.7];    P(2,6,:)=[1 3 0];   P(2,7,:)=[1 4 0.8];
P(3,1,:)=[2 0 0.8];   P(3,2,:)=[2 1 0];  P(3,3,:)=[2 2 0.7];P(3,4,:)=[2 2 0.7];P(3,5,:)=[2 2 0.7];      P(3,6,:)=[2 3 0];   P(3,7,:)=[2 4 0.8];
P(4,1,:)=[3 0 0.6];   P(4,2,:)=[3 1 0];   P(4,3,:)=[3 2 0.5];P(4,4,:)=[3 2 0.5];P(4,5,:)=[3 2 0.5];     P(4,6,:)=[3 3 0];   P(4,7,:)=[3 4 0.5];
P(5,1,:)=[4 0 0.8];   P(5,2,:)=[4 1 0];   P(5,3,:)=[4 2 0.7];P(5,4,:)=[4 2 0.7];P(5,5,:)=[4 2 0.7];     P(5,6,:)=[4 3 0];   P(5,7,:)=[4 4 0.8];
P(6,1,:)=[5 0 0.8];   P(6,2,:)=[5 1 0];    P(6,3,:)=[5 2 0.7]; P(6,4,:)=[5 2 0.7];P(6,5,:)=[5 2 0.7];   P(6,6,:)=[5 3 0];   P(6,7,:)=[5 4 0.8];
P(7,1,:)=[6 0 1];   P(7,2,:)=[6 1 1];    P(7,3,:)=[6 2 0.8];  P(7,4,:)=[6 2 0.8];P(7,5,:)=[6 2 0.8];          P(7,6,:)=[6 3 1];     P(7,7,:)=[6 4 1];
figure
plot3(P(1:7,:,1),P(1:7,:,2),P(:,:,3));xlabel('x');ylabel('y');grid on;
figure
surf(P(1:7,:,1),P(1:7,:,2),P(:,:,3),'FaceAlpha',0.5)
% i=numero elementi su x; j n ele su y

%% automatic cefinition of equispaced knot 
pu=2; % degree of curve u
U=zeros(1,numel(P(1,:,1))+pu+1);
U(end-pu:end)=1;
mid=numel(U)-2*pu-2;
for i=1:mid
    U(i+pu+1)=1/(mid+1)*i;
end

pv=2;
V=zeros(1,numel(P(:,1,1))+pv+1);
V(end-pv:end)=1;
mid=numel(V)-2*pv-2;
for i=1:mid
    V(i+pv+1)=1/(mid+1)*i;
end

u=0:0.01:1-0.01; v=u;
%% surface creation
%for iu=1:numel(u)
u=0.5;
    for iv=1:numel(v)
        S(iv)=0;
        for i=1:numel(P(:,1,2))
            for j=1:numel(P(1,:,1))
                S(iv)=S(iv)+Nip(j,u,pu,U)*Nip(i,v(iv),pv,V)*P(i,j,3);
            end
        end
    end
%end
figure 
rx=max(max(P(:,:,1))); ry=max(max(P(:,:,2)));
plot(v*ry,S);hold on
surf(P(1:5,:,1),P(:,1:5,2),P(:,:,3),'FaceAlpha',0,'MarkerFaceColor','r','LineStyle','none','Marker','o','MarkerSize',7)
