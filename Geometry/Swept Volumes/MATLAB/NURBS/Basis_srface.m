clear all
close all
clc
% i span, nel libro i parte da 0, qua da 1
U=[0 0 0 1 2 3 4 5 5 5];
V=[0 0 1 1];
lowv=2; highv=numel(V)-lowv-1;
lowu=2; highu=numel(U)-lowu-1;

u=min(U):0.1:max(U);
v=u;
F=meshgrid(u,v);
for pv=0:lowv-1
    for iv=1:highv
        for pu=0:lowu-1
                for iu=1:highu        
                    for j=1:numel(v)
                        Nipv=Nip(iv,v(j),pv,V);
                        for k=1:numel(u) 
                            N(j,k)=Nip(iu,u(k),pu,U)*Nipv;
                        end
                    end
                    figure
                    surf(u,v,N);title(['iv=',num2str(iv),' ; pv=',num2str(pv),' iu=',num2str(iu),' ; pu=',num2str(pu)])
            end
        end
        
    end
   
end
%% surface from point
clear all
P(1,:,:)=[3 0 1 ; 1 0 1; 2 0 3; 3 0 3]';
P(2,:,:)=[3 1 1 ; 1 1 0; 2 1 2; 3 1 2]';
P(3,:,:)=[3 2 1 ; 1 2 0; 2 2 2; 3 2 2]';

n=4;  m=3;

U=[0 0 0 1 4 5 5 5];
V=[0 0 0 2 3 3 3];% more molteplicity more continuity :)
lowv=2; highv=numel(V)-lowv-1;
lowu=2; highu=numel(U)-lowu-1;
u=0:0.1:5-0.1;
v=0:0.1:3-0.1;
p=2;

for j=1:numel(u)
    for k=1:numel(v)
        S(j,k)=0;
        for i=1:n
            for l=1:m
                S(j,k)=S(j,k)+Nip(i,u(j),p,U)*Nip(l,v(k),p,V)*P(l,3,i); 
            end
        end
    end
end
surf(v,u,S)
    



