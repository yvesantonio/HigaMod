% volCoefs = zeros(4,N,N,N);
% clear all
% close all
% clc
% addpath('C:\Users\XPS\Desktop\PC port\NURBS',...
% 'C:\Users\XPS\Desktop\PC port\nurbs_toolbox',...
% 'C:\Users\XPS\Desktop\PC port\My_Function');
% 

% load('section_model_1_1_380.mat')
%% get coons surface in nurbs structure
for i=1:numel(surf_point)
    coons{i}=get_coons(surf_point{i},10,2);
end
%%
k=1;
for i = 1:numel(coons)
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
nrbplot(VOL,[10 10 150]);
