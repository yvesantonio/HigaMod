clear all
close all
clc
%% definizione superficie tramite punti versp il basso X; verso destra Y; P(x,y,z)
P(1,1,:)=[0 0 0];   P(1,2,:)=[1 0 0];   P(1,3,:)=[2 0 0];   P(1,4,:)=[3 0 0];    P(1,5,:)=[4 0 0];   P(1,6,:)=[5 0 0];   P(1,7,:)=[6 0 0];
P(2,1,:)=[0 1 0];  P(2,2,:)=[1 1 1];   P(2,3,:)=[2 1 0];   P(2,4,:)=[3 1 1];    P(2,5,:)=[4 1 0];   P(2,6,:)=[5 1 0];   P(2,7,:)=[6 1 0];
P(3,1,:)=[0 2 0];   P(3,2,:)=[1 2 0];   P(3,3,:)=[2 2 1];  P(3,4,:)=[3 2 0];     P(3,5,:)=[4 2 0];  P(3,6,:)=[5 2 0.0]; P(3,7,:)=[6 2 0.0];
P(4,1,:)=[0 3 1];   P(4,2,:)=[1 3 0];   P(4,3,:)=[2 3 0];   P(4,4,:)=[3 3 1];    P(4,5,:)=[4 3 0]; P(4,6,:)=[5 3 0.0]; P(4,7,:)=[6 3 0.0];
P(5,1,:)=[0 4 1];   P(5,2,:)=[1 4 0];   P(5,3,:)=[2 4 0];   P(5,4,:)=[3 4 1];  P(5,5,:)=[4 4 0];  P(5,6,:)=[5 4 0];   P(5,7,:)=[6 4 0];
P(6,1,:)=[0 5 0];   P(6,2,:)=[1 5 1];   P(6,3,:)=[2 5 1];   P(6,4,:)=[3 5 -0];  P(6,5,:)=[4 5 0];  P(6,6,:)=[5 5 0];   P(6,7,:)=[6 5 0];
P(7,1,:)=[0 6 0];   P(7,2,:)=[1 6 0];   P(7,3,:)=[2 6 1];   P(7,4,:)=[3 6 -0];  P(7,5,:)=[4 6 0];  P(7,6,:)=[5 6 0];   P(7,7,:)=[6 6 0];
P(8,1,:)=[0 7 0];   P(8,2,:)=[1 7 0];   P(8,3,:)=[2 7 0];   P(8,4,:)=[3 7 -0];  P(8,5,:)=[4 7 0];  P(8,6,:)=[5 7 0];   P(8,7,:)=[6 7 0];
P(9,1,:)=[0 8 0];   P(9,2,:)=[1 8 0];   P(9,3,:)=[2 8 0];   P(9,4,:)=[3 8 -0];  P(9,5,:)=[4 8 0];   P(9,6,:)=[5 8 0];  P(9,7,:)=[6 8 0];
%% riscalo x=7, y=9

W=[1 1 1 1 1 1 1;
   1 1 1 1 1 1 1;
   1 1 1 1 1  1 1;
   1 1 1 1 1 1 1;
   1 1 1 1 1 1 1;
    1 1 1 1 1 1 1;
     1 1 1 1 1 1 1;
      1 1 1 1 1 1 1;
      1 1 1 1 1 1 1];
figure
plot3(P(1:9,:,1),P(:,1:7,2),P(:,:,3),'rd');xlabel('x');ylabel('y');grid on
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

u=[0:0.01:1-0.01]; v=u;
%% surface creation
for iu=1:numel(u)
    for iv=1:numel(v)
        S(iv,iu,:)=[0;0;0];
        for i=1:numel(P(:,1,2))
            for j=1:numel(P(1,:,1))
                S(iv,iu,:)=S(iv,iu,:)+Nip(j,u(iu),pu,U)*Nip(i,v(iv),pv,V)*P(i,j,:)*W(i,j);
            end
        end
        div=0;
        for i=1:numel(P(:,1,2))
            for j=1:numel(P(1,:,1))
                div=div+Nip(j,u(iu),pu,U)*Nip(i,v(iv),pv,V)*W(i,j);
            end
        end
       S(iv,iu,:)=S(iv,iu,:)/div;
    end
end
figure 
axis('image');
rx=max(max(P(:,:,1))); ry=max(max(P(:,:,2)));
surf(u*rx,v*ry,S(:,:,3),'FaceAlpha',1);hold on
surf(P(1,:,1),P(:,1,2),P(:,:,3),'FaceAlpha',0,'MarkerFaceColor','r','LineStyle','none','Marker','o','MarkerSize',7)
figure
plot3(   S(1:100,1:100,1), S(1:100,1:100,2),S(1:100,1:100,3),'b'  );hold on
