clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');


load('point_dist_prova.mat')
load('point_coordtg.mat')
load('point_provatg.mat')
%% mettere tutto in senso di percorrenza orario
for i=1:numel(surf_point)
ss=surf_point{i};
vett1=surf_point{i}(40,:)-coordtg(i,:);
vett2=surf_point{i}(80,:)-coordtg(i,:);
veri=[1;0;0]; verj=[0;1;0];  verk=[0;0;1];
prodvett=veri*det([vett1(2) vett1(3);vett2(2) vett2(3)])-...
         verj*det([vett1(1) vett1(3);vett2(1) vett2(3)])+...
         verk*det([vett1(1) vett1(2);vett2(1) vett2(2)]);
tg=Tgsrf(i,:);
prodscal(i)=tg*prodvett;
if prodscal(i)<0  %inverti il senso antior con or
   for j=1:numel(ss)/3
       ssnew(j,:)=ss(end-j+1,:);
   end
   surf_point{i}=ssnew;
end

end
min(prodscal)
%%
for i=1:numel(surf_point)
ss=surf_point{i};
vett1=surf_point{i}(40,:)-coordtg(i,:);
vett2=surf_point{i}(80,:)-coordtg(i,:);
veri=[1;0;0]; verj=[0;1;0];  verk=[0;0;1];
prodvett=veri*det([vett1(2) vett1(3);vett2(2) vett2(3)])-...
         verj*det([vett1(1) vett1(3);vett2(1) vett2(3)])+...
         verk*det([vett1(1) vett1(2);vett2(1) vett2(2)]);
tg=Tgsrf(i,:);
prodscal(i)=tg*prodvett;
% if prodscal<0  %inverti il senso antior con or
%    for j=1:numel(ss)/3
%        ssnew(j,:)=ss(end-j+1,:);
%    end
% end
end
min(prodscal)
%%
figure
hold on
grid on
srfn{1}=surf_point{1};
for k=2:numel(surf_point)
    
srf=surf_point{k};
nel=numel(srf)/3;
df=norm(srf(1,:));
dn=norm(srf(1,:));
far=1;
near=1;
for i=2:nel
     if norm(srf(i,3))>df
         df=norm(srf(i,:));
         far=i;
     end
     if norm(srf(i,:))<dn
         dn=norm(srf(i,:));
         near=i;
     end
end 
clearvars xn yn zn
xn=zeros(nel,1); yn=zeros(nel,1); zn=zeros(nel,1);
if far-near < 0
    xn(1:nel-far+1,1)=srf(far:end,1); yn(1:nel-far+1,1)=srf(far:end,2); zn(1:nel-far+1,1)=srf(far:end,3);
    xn(nel-far+2:end,1)=srf(2:far,1);  yn(nel-far+2:end,1)=srf(2:far,2); zn(nel-far+2:end,1)=srf(2:far,3); 
else
    xn(1:far)=srf(far:-1:1,1); yn(1:far)=srf(far:-1:1,2); zn(1:far)=srf(far:-1:1,3);
    xn(far+1:end,1)=srf(end-1:-1:far,1); yn(far+1:end,1)=srf(end-1:-1:far,2); zn(far+1:end,1)=srf(end-1:-1:far,3); 
end
srfn{k}=[xn yn zn];
% plot3(srf(far,1),srf(far,2),srf(far,3),'rd',srf(far+10,1),srf(far+10,2),srf(far+10,3),'gd')
plot3(xn,yn,zn,'b',xn(1),yn(1),zn(1),'rd',xn(10),yn(10),zn(10),'gd')  %,srf(near,1),srf(near,2),srf(near,3),'bd'
end

%% get coons surface in nurbs structure
for i=1:numel(srfn)
    coons{i}=get_coons(srfn{i},6,2);
end
%%
k=1;
for i = 1:numel(srfn)
    surf_i = coons{i};
    volCoefs(:,:,:,k) = surf_i.coefs;
    k=k+1;
end
%% deg of volume 
d=2;
%% on U direction
n=numel(volCoefs(1,:,1,1));
U=zeros(1,n+d+1); % automatic node creation
U(end-d:end)=1;
mid=numel(U)-2*d-2;
for i=1:mid
    U(i+d+1)=1/(mid+1)*i;
end
%% on V direction
n=numel(volCoefs(1,1,:,1));
V=zeros(1,n+d+1); % automatic node creation
V(end-d:end)=1;
mid=numel(V)-2*d-2;
for i=1:mid
    V(i+d+1)=1/(mid+1)*i;
end
%% on W direction
n=numel(volCoefs(1,1,1,:));
W=zeros(1,n+d+1); % automatic node creation
W(end-d:end)=1;
mid=numel(W)-2*d-2;
for i=1:mid
    W(i+d+1)=1/(mid+1)*i;
end
%%
VOL = nrbmak(volCoefs,{U,V,W});
figure
nrbplot(VOL,[10 10 50]);
