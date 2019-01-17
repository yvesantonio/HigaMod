function [ best_approx_eval, best_approx  ]=g_n_opt(u,x1,y1,z1,x,y,z,U,tol,max_k,max_k_l)
n=numel(x1);
for i=1:numel(x1)
    pL(i*3-2:i*3,1)=[x1(i);y1(i);z1(i)];
end

x_LS=[u';pL;ones(n,1)]; % starting guess for optimal solution
m=numel(u);
k=1;
% per entrrare nel ciclo per forza
norm_temp=inf;
norm_fu=0;
while k<=max_k & abs(norm_temp-norm_fu)>tol
    
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
    Jaug=[J; zeros(n,m+3*n) a*eye(n)];  %
    fuaug=[fu;a*(x_LS(end-n+1:end)-1)]; %
    % solution to find p
    p=mygaus((Jaug'*Jaug),(-Jaug'*fuaug),zeros(m+4*n,1),1.e-10,1000);

    % set lambda  l = 1 and start "line search"
    l=1;
    i_l=1;
    go_on=1;
    while go_on & i_l<=max_k_l
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
    k=k+1;
hold on
% debug verificare convergenza con disegno e numeri
% plot3(best_approx_eval(1,:),best_approx_eval(2,:),best_approx_eval(3,:),'Color',[1-lim_k/20 0 lim_k/20])
% plot3(x_(m+1:3:m+3*n)',x_(m+2:3:m+3*n)',x_(m+3:3:m+3*n),'d','Color',[1-lim_k/20 0 lim_k/20])
%         norm_temp=norm(fu_temp)
%         norm_fu=norm(fu)                                            
end
% forzare inizio e fine a rispettare la curva
% x_LS(m+1:m+3,1)=[x(1);y(1);z(1)];
% x_LS(end-n-2:end-n)=[x(end);y(end);z(end)];
best_approx=nrbmak([x_LS(m+1:3:m+3*n)';...
                 x_LS(m+2:3:m+3*n)';...
                 x_LS(m+3:3:m+3*n)';...
                 x_LS(end-n+1:end)'],U);
best_approx_eval=nrbeval(best_approx,x_LS(1:m));
% best_approx_eval(1:3,1)=[x(1);y(1);z(1)];
% best_approx_eval(1:3,end)=[x(end);y(end);z(end)];
end