x_LS=[u';pL;ones(n,1)]; % da inizializzazione
m=numel(u);
k=1;
lim_k=5;
norm_temp=inf;
norm_fu=0;
while k<=lim_k & abs(norm_temp-norm_fu)>1.e-4
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
    Jaug=[J; zeros(n,m+3*n) a*eye(n)];  %
    fuaug=[fu;a*(x_LS(end-n+1:end)-1)]; %
    % solution to find p
    p=mygaus((Jaug'*Jaug),(-Jaug'*fuaug),zeros(m+4*n,1),1.e-5,1000);

    % set lambda  l = 1 and start "line search"
    l=1;
    i_l=1;
    go_on=1;
    while go_on & i_l<=20
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
plot3(best_approx_eval(1,:),best_approx_eval(2,:),best_approx_eval(3,:),'Color',[1-lim_k/20 0 lim_k/20])
plot3(x_(m+1:3:m+3*n)',x_(m+2:3:m+3*n)',x_(m+3:3:m+3*n),'d','Color',[1-lim_k/20 0 lim_k/20])
        norm_temp=norm(fu_temp)
        norm_fu=norm(fu)                                            
end

%% forzare gli ultimi punti ad essere i punti
x_LS(m+1:m+3,1)=[x(1);y(1);z(1)];
x_LS(end-n-2:end-n)=[x(end);y(end);z(end)];
best_approx=nrbmak([x_LS(m+1:3:m+3*n)';...
                 x_LS(m+2:3:m+3*n)';...
                 x_LS(m+3:3:m+3*n)';...
                 x_LS(end-n+1:end)'],U);
best_approx_eval=nrbeval(best_approx,x_LS(1:m));
best_approx_eval(end-2:end)=[x(end);y(end);z(end)];

%% conronto finale perchè ne so
figure; hold on
plot3(x,y,z,'r'); grid on; title('curva iniziale in r, firt approx in b, finish in k')
best_approx=nrbmak([x_LS(m+1:3:m+3*n)';...
                 x_LS(m+2:3:m+3*n)';...
                 x_LS(m+3:3:m+3*n)';...
                 x_LS(end-n+1:end)'],U);
best_approx_eval=nrbeval(best_approx,x_LS(1:m));

best_approx_eval(end-2:end)=[x(end);y(end);z(end)];
plot3(best_approx_eval(1,:),best_approx_eval(2,:),best_approx_eval(3,:),'k')
plot3(x_LS(m+1:3:m+3*n)',x_LS(m+2:3:m+3*n)',x_LS(m+3:3:m+3*n),'kd')
first_approx=nrbmak([x1';y1';z1'],U);
nrbplot(first_approx,100)
