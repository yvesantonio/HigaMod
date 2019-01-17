%% first iteration to find the starting control point (from other script)

sol_LS=[u';pL;w]; % simple LS with all w=1
k=1;
while k==1 | k<=5
    k
    % define f(u)
    nrb_now=nrbmak([sol_LS(numel(u)+1:3:numel(u)+3*n)';sol_LS(numel(u)+2:3:numel(u)+3*n)';...
        sol_LS(numel(u)+3:3:numel(u)+3*n)';sol_LS(end-n+1:end)'],U);
        P = nrbeval(nrb_now,[0:0.05:0.95 0.9999]);
    for i=1:numel(u)
        fu(i*3-2:i*3,1)=[P(1,i)-x(i);P(2,i)-y(i);P(3,i)-z(i)];
    end
%     a=norm(fu)^2/n;
%     % computation of J_aug and f_aug
    J=get_J(nrb_now,u,n);
%     Jaug=[J; zeros(n,m+3*n) a*eye(n)];  %
%     fuaug=[fu;a*(sol_LS(end-n+1:end)-1)]; %
    % solution to find p
    p=(J'*J)\(J'*fu);
    % set lambda  l = 1 and start "line search"
    l=1;
    go_on=1;
    while go_on
        pl=p*l; 
        x_=sol_LS+pl;
        x_=bound_apply(x_,u,w);
        % valutazione della norma dopo il passo
        nrb_temp=nrbmak([x_(numel(u)+1:3:numel(u)+3*n)';x_(numel(u)+2:3:numel(u)+3*n)';...
                        x_(numel(u)+3:3:numel(u)+3*n)';x_(end-n+1:end)'],U);
        nrb_tempeval=nrbeval(nrb_temp,x_(1:numel(u)));
        for i=1:numel(u)
            fu_temp(i*3-2:i*3,1)=[nrb_tempeval(1,i)-x(i);nrb_tempeval(2,i)-y(i);nrb_tempeval(3,i)-z(i)];
        end
        norm_temp=norm(fu_temp)
        norm_fu=norm(fu)
        if  norm_temp < norm_fu
            go_on=0;
            sol_LS=x_;
        end
        l=l/2;
        i_l=i_l+1;
    end
    sol_LS=x_;
    k=k+1;
end
%% confronto grafico
hold on
nrb_now=nrbmak([sol_LS(numel(u)+1:3:numel(u)+3*n)';sol_LS(numel(u)+2:3:numel(u)+3*n)';...
        sol_LS(numel(u)+3:3:numel(u)+3*n)';sol_LS(end-n+1:end)'],U);
        P = nrbeval(nrb_now,[0:0.05:0.95 0.9999]);
plot3(P(1,:),P(2,:),P(3,:),'--g')
%% confronto numerico
