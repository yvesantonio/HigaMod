clear all
close all
clc
addpath('C:\Users\XPS\Desktop\PC port\NURBS',...
'C:\Users\XPS\Desktop\PC port\nurbs_toolbox',...
'C:\Users\XPS\Desktop\PC port\My_Function');

load('model_to_section_yves.mat')
%% " rotate " to minimiza the distance
% s1=surf_point{1};s2=surf_point{2};
% 
% [x y z]=nearest(s1(1,1),s1(1,2),s1(1,3),s2(:,1),s2(:,2),s2(:,3));
% [x2 y2 z2]=nearest(s1(2,1),s1(2,2),s1(2,3),s2(:,1),s2(:,2),s2(:,3));
%% tutte le surf sone definite dallo stesso numero di punti
figure
hold on
grid on
npti=numel(surf_point{1}(:,1));
for i=40:50
    npti=max(npti,numel(surf_point{i}(:,1)));
end
ev=linspace(0,1,npti);
for i=40:50
    plot3(surf_point{i}(:,1),surf_point{i}(:,2),surf_point{i}(:,3),'r')
    [a nrb_temp]=nrb_approx(surf_point{i}(:,1),surf_point{i}(:,2),surf_point{i}(:,3),ceil(npti/10),3,1.e-5,20,100);
    new_srf{i-39}=nrbeval(nrb_temp,ev)';
    new_srf{i-39}(end,:)=new_srf{i-39}(1,:);
    plot3(new_srf{i-39}(:,1),new_srf{i-39}(:,2),new_srf{i-39}(:,3),'b')
    pause(1)
end

% %% riordina i punti... -- da estendere a funzione -- 
% x1=new_srf{1}(:,1); y1=new_srf{1}(:,2); z1=new_srf{1}(:,3);
% x2=new_srf{2}(:,1); y2=new_srf{2}(:,2); z2=new_srf{2}(:,3);
% xn=zeros(numel(x2),1); yn=zeros(numel(x2),1); zn=zeros(numel(x2),1); 
% [xn(1) yn(1) zn(1)]=nearest(x1(1),y1(1),z1(1),x2,y2,z2);
% xn(end)=xn(1); yn(end)=yn(1); zn(end)=zn(1); 
% for i=2:npti-1
%     xn(i)=inf; yn(i)=inf; zn(i)=inf; 
%     for j=1:npti-1
%         isnot=isnotin(x2(j),y2(j),z2(j),xn,yn,zn);
%         if isnot & norm([x1(i) y1(i) z1(i)]-[xn(i) yn(i) zn(i)])>norm([x1(i) y1(i) z1(i)]-[x2(j) y2(j) z2(j)])
%                    
%                xn(i)=x2(j); yn(i)=y2(j); zn(i)=z2(j); 
%         end
%     end
% 
%                
% end
%%
figure
hold on
grid on
for i=1:numel(new_srf)
    coons_srf(i)=get_coons(new_srf{i},10,2);
   
    nrbplot(coons_srf(i),[30 30])
end
          














