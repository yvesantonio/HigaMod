%% initialize xLS vector of the solution
xLS=[u';pL;w]; % simple LS with all w=1
nrb=nrbmak([xL';yL';zL'],U);
P=nrbeval(nrb,[0:0.05:0.95 0.9999]);
% create fu vector
for i=1:numel(u)
fu(i*3-2:i*3,1)=[P(1,i)-x(i);P(2,i)-y(i);P(3,i)-z(i)];
end
%% initialize G-N with line search
tol=1.e-4;
i=1; max_iter=20;
a=norm(fu)^2/n;
J=get_J(nrb,u,n);
Jaug=[J; zeros(n,m+3*n) a*eye(n)];
fuaug=[fu;a*xLS(end-n+1:end)];
p=(Jaug'*Jaug)\(-Jaug'*fuaug);
% line search
l=1;
x_=xLS+p*l;
x_=bound_apply(x_,u,w);
l=l/2;
nrb_temp=nrbmak([x_(numel(u)+1:3:numel(u)+3*n)';x_(numel(u)+2:3:numel(u)+3*n)';x_(numel(u)+3:3:numel(u)+3*n)';x_(end-n+1:end)'],U);
nrb_tempe=nrbeval(nrb_temp,x_(1:numel(u)));
for i=1:numel(u)
fu_temp(i*3-2:i*3,1)=[nrb_tempe(1,i)-x(i);nrb_tempe(2,i)-y(i);nrb_tempe(3,i)-z(i)];
end
norm(fu_temp)
norm(fu)
while norm(fu_temp)>=norm(fu);
    x_=xLS+p*l;
    x_=bound_apply(x_,u,w);
    l=l/2;
    nrb_temp=nrbmak([x_(numel(u)+1:3:numel(u)+3*n)';x_(numel(u)+2:3:numel(u)+3*n)';x_(numel(u)+3:3:numel(u)+3*n)';x_(end-n+1:end)'],U);
    nrb_tempe=nrbeval(nrb_temp,x_(1:numel(u)));
    for i=1:numel(u)
    fu_temp(i*3-2:i*3,1)=[nrb_tempe(1,i)-x(i);nrb_tempe(2,i)-y(i);nrb_tempe(3,i)-z(i)];
    end
    norm(fu_temp)
    norm(fu)
end
xLS=x_;
while i<=max_iter
    nrb=nrbmak([xLS(numel(u)+1:3:numel(u)+3*n)';xLS(numel(u)+2:3:numel(u)+3*n)';...
        xLS(numel(u)+3:3:numel(u)+3*n)';xLS(end-n+1:end)'],U);
    for i=1:numel(u)
        fu(i*3-2:i*3,1)=[P(1,i)-x(i);P(2,i)-y(i);P(3,i)-z(i)];
    end

    
    
end
    
    
    
    
    