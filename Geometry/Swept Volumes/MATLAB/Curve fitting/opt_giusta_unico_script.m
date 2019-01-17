clear all
close all
clc
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
  'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
  'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%%
% x=[0:0.01:2];
% y=[sin(x*3*pi)+sqrt(x)];
% z=zeros(1,numel(x));
% figure; plot3(x,y,z,'r'); grid on; title('curva iniziale')
x=[0:0.01:2 2.01:0.01:3.99 4:0.01:6];
z=zeros(1,numel(x));
y=[cos(4*(0:0.01:2).^2) -sin(exp((0:0.01:1.98))) 4-2*sin(((4:0.01:6)-4))];
figure; plot3(x,y,z); grid on; title('curva iniziale')
%% fitting per i punti di partenza con pesi = 1
n=35;
d=3;
u=linspace(0,1,numel(x)); u(end)=0.99999999999;
m=numel(u);
U=zeros(1,n+d+1); % automatic node creation
mid=numel(U)-2*d-2;
for i=1:2:mid
    U(i+d+1:i+d+2)=[1/(mid+1)*i ];
end
U(end-d:end)=1;
%% u(1)=0 fisso u(m)=0.99999 fisso
clearvars element
for i=1:n
    for j=1:m
    element(j,i)=Nip(i,u(j),d,U);
    end
end


bx=x';
by=y';
bz=z';
x1=element\bx;
y1=element\by;
z1=element\bz;

% x1(1)=x(1); x1(n)=x(end);
% y1(1)=y(1); y1(n)=y(end);
% z1(1)=z(1); z1(n)=z(end);
hold on; plot3(x1,y1,z1,'bd')
%%
first_approx=nrbmak([x1';y1';z1'],U);
nrbplot(first_approx,100)
for i=1:numel(x1)
    pL(i*3-2:i*3,1)=[x1(i);y1(i);z1(i)];
end
%% Ju for gauss-N
fa_der=nrbderiv(first_approx);
 [pnt, jac] = nrbdeval (first_approx, fa_der, u);
Ju=zeros(3*numel(u),numel(u));
for i=1:numel(u)
    Ju(i*3-2:i*3,i)=jac(:,i);
end
%% Jp for gauss_n
Jp=zeros(3*numel(u),3*n);
for i=1:numel(u)
    for j=1:first_approx.number;
    der(i,j)=nrbeval_der_p(first_approx, j,u(i));
    Jp(i*3-2:i*3,j*3-2:j*3)=diag([der(i,j),der(i,j),der(i,j)]);
    end
 end
%% Jw fo G-N
Jw=zeros(3*numel(u),n);
for i=1:numel(u)
    for j=1:first_approx.number;
    Jw(i*3-2:i*3,j)=nrbeval_der_w(first_approx, j,u(i)) ;
    end
end
%% create J
J=[Ju Jp Jw];

%%
figure; hold on
x_LS=[u';pL;ones(n,1)]; % da inizializzazione
m=numel(u);
k=1;
lim_k=20;
norm_temp=inf;
norm_fu=0;
while k<=lim_k & abs(norm_temp-norm_fu)>1.e-10
    k % vedo dove caxxo sono
    % define f(u)
    nrb_now=nrbmak([x_LS(m+1:3:m+3*n)';... % struttura nurbs con quello che penso 
                    x_LS(m+2:3:m+3*n)';... % sia la soluzione migliore
                    x_LS(m+3:3:m+3*n)';...
                    x_LS(end-n+1:end)'],U); 
    P = nrbeval(nrb_now,x_LS(1:m));    % valutazione della nurbs nei punti
    for i=1:numel(u)
        fu(i*3-2:i*3,1)=[P(1,i)-x(i);P(2,i)-y(i);P(3,i)-z(i)];
    end  % fu() --> differenza tra punti ora e punti misurati
    a=norm(fu)^2/k;
    
    % computation of J_aug and f_aug
    
    J=get_J(nrb_now,x_LS(1:m),n);
     Jaug=[J ;zeros(n,m+3*n) a*eye(n)];  % 
     fuaug=[fu;a*(x_LS(end-n+1:end)-1)]; % 
    % solution to find p
    p=mygaus((Jaug'*Jaug),(-Jaug'*fuaug),zeros(m+4*n,1),1.e-5,1000);

    % set lambda  l = 1 and start "line search"
    l=1;
    i_l=1;
    go_on=1;
    while go_on & i_l<=1000
        pl=p*l; 
        x_=x_LS+pl;
        x_=bound_apply(x_,m,n);
        % valutazione della norma dopo il passo
        best_approx=nrbmak([x_(m+1:3:m+3*n)';...
                         x_(m+2:3:m+3*n)';...
                         x_(m+3:3:m+3*n)';...
                         x_(end-n+1:end)'],U);
        best_approx_eval=nrbeval(best_approx,x_(1:m));
        for i=1:numel(u)
            fu_temp(i*3-2:i*3,1)=[best_approx_eval(1,i)-x(i);best_approx_eval(2,i)-y(i);best_approx_eval(3,i)-z(i)];
        end
        norm_temp=norm(fu_temp);
        norm_fu=norm(fu);
        if  norm_temp < norm_fu
            go_on=0;
        end
        l=l/2;
        i_l=i_l+1;
    end
    x_LS=x_;
   
hold on
plot3(best_approx_eval(1,:),best_approx_eval(2,:),best_approx_eval(3,:),'Color',[1 0 1],'LineWidth',0.5)
% plot3(x_(m+1:3:m+3*n)',x_(m+2:3:m+3*n)',x_(m+3:3:m+3*n),'d','Color',[1-lim_k/20 0 lim_k/20])
        norm_temp=norm(fu_temp)
        norm_fu=norm(fu) 
        delta(k)=abs(norm_temp-norm_fu);
         k=k+1;
end

%% forzare gli ultimi punti ad essere i punti
%x_LS(m+1:m+3,1)=[x(1);y(1);z(1)];
%x_LS(end-n-2:end-n)=[x(end);y(end);z(end)];
best_approx=nrbmak([x_LS(m+1:3:m+3*n)';...
                 x_LS(m+2:3:m+3*n)';...
                 x_LS(m+3:3:m+3*n)';...
                 x_LS(end-n+1:end)'],U);
%best_approx_eval=nrbeval(best_approx,x_LS(1:m));
%best_approx_eval(end-2:end)=[x(end);y(end);z(end)];

%% conronto finale perchè ne so
figure; hold on
plot3(x,y,z,'r'); grid on; title('curva iniziale in r, firt approx in b, finish in k')
% best_approx=nrbmak([x_LS(m+1:3:m+3*n)';...
%                  x_LS(m+2:3:m+3*n)';...
%                  x_LS(m+3:3:m+3*n)';...
%                  x_LS(end-n+1:end)'],U);
% best_approx_eval=nrbeval(best_approx,x_LS(1:m));

%best_approx_eval(end-2:end)=[x(end);y(end);z(end)];
plot3(best_approx_eval(1,:),best_approx_eval(2,:),best_approx_eval(3,:),'k')
plot3(x_LS(m+1:3:m+3*n)',x_LS(m+2:3:m+3*n)',x_LS(m+3:3:m+3*n),'kd')
first_approx=nrbmak([x1';y1';z1'],U);
plot3(x1,y1,z1,'rd')
nrbplot(first_approx,100)
delta


