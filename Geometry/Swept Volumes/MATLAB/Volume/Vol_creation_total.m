clear all
close all
clc
% addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
% 'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
% 'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');


load('coord_total.mat')
load('section_model_total.mat')
load('Tg_total.mat')
%% mettere tutto in senso di percorrenza orario
for i=1:numel(surf_point)
ss=surf_point{i};
vett1=surf_point{i}(1,:)-coordtg(i,:);
vett2=surf_point{i}(round(numel(ss)/3/4),:)-coordtg(i,:);
veri=[1;0;0]; verj=[0;1;0];  verk=[0;0;1];
prodvett=veri*det([vett1(2) vett1(3);vett2(2) vett2(3)])-...
         verj*det([vett1(1) vett1(3);vett2(1) vett2(3)])+...
         verk*det([vett1(1) vett1(2);vett2(1) vett2(2)]);
tg=tgsrf(i,:);
prodscal_best(i)=tg*prodvett;
if prodscal_best(i)<0  %inverti il senso antior con or
   for j=1:numel(ss)/3
       ssnew(j,:)=ss(end-j+1,:);
   end
   surf_point{i}=ssnew;
   clearvars ssnew
end

end
min(prodscal_best)
%%
for i=1:numel(surf_point)
ss=surf_point{i};
vett1=surf_point{i}(1,:)-coordtg(i,:);
vett2=surf_point{i}(20,:)-coordtg(i,:);
veri=[1;0;0]; verj=[0;1;0];  verk=[0;0;1];
prodvett=veri*det([vett1(2) vett1(3);vett2(2) vett2(3)])-...
         verj*det([vett1(1) vett1(3);vett2(1) vett2(3)])+...
         verk*det([vett1(1) vett1(2);vett2(1) vett2(2)]);
tg=tgsrf(i,:);
prodscal_best(i)=tg*prodvett;
% if prodscal<0  %inverti il senso antior con or
%    for j=1:numel(ss)/3
%        ssnew(j,:)=ss(end-j+1,:);
%    end
% end
end
min(prodscal_best)
%%
figure
hold on
grid on
for i=1:numel(surf_point)

plot3(surf_point{i}(:,1),surf_point{i}(:,2),surf_point{i}(:,3),'b',...
      surf_point{i}(1,1),surf_point{i}(1,2),surf_point{i}(1,3),'rd',...
      surf_point{i}(10,1),surf_point{i}(10,2),surf_point{i}(10,3),'gd')
end
%% oridane
figure
hold on
grid on
srfn{1}=surf_point{1};
for k=2:numel(surf_point)
    srf=surf_point{k};
    vett_prec=srfn{k-1}(1,:)-coordtg(k-1,:); vett_prec=vett_prec/norm(vett_prec);
    nel=numel(surf_point{k})/3;
    best=1;
    worst=1;
    prodscal_best=vett_prec*((surf_point{k}(1,:)-coordtg(k,:))'/norm((surf_point{k}(1,:)-coordtg(k,:))'));
    prodscal_worst=vett_prec*((surf_point{k}(1,:)-coordtg(k,:))'/norm((surf_point{k}(1,:)-coordtg(k,:))'));
    for j=2:nel-1
         guess=(surf_point{k}(j,:)-coordtg(k,:))/norm(surf_point{k}(j,:)-coordtg(k,:));
         prodscal_new=vett_prec*guess';
         if prodscal_new > prodscal_best
             best=j;
             prodscal_best=prodscal_new;
         end
         if prodscal_new < prodscal_worst
             worst=j;
             prodscal_worst=prodscal_new;
         end
    end
    clearvars xn yn zn
    xn=zeros(nel,1); yn=zeros(nel,1); zn=zeros(nel,1);
        xn(1:nel-best+1,1)=srf(best:end,1); yn(1:nel-best+1,1)=srf(best:end,2); zn(1:nel-best+1,1)=srf(best:end,3);
        xn(nel-best+2:end,1)=srf(2:best,1);  yn(nel-best+2:end,1)=srf(2:best,2); zn(nel-best+2:end,1)=srf(2:best,3);  

    srfn{k}=[xn yn zn];
    % plot3(srf(far,1),srf(far,2),srf(far,3),'rd',srf(far+10,1),srf(far+10,2),srf(far+10,3),'gd')
    plot3(xn,yn,zn,'k',xn(1),yn(1),zn(1),'kd')  %,srf(near,1),srf(near,2),srf(near,3),'bd'
end
%% get coons surface in nurbs structure
for i=1:numel(srfn)
    coons{i}=get_coons(srfn{i},4,2);
end
%%
k=1;
for i = 1:numel(srfn)
    surf_i = coons{i};
    volCoefs(:,:,:,k) = surf_i.coefs;
    k=k+1;
end
%% deg of volume 
d=3;
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
nrbplot(VOL,[20 20 800]);
