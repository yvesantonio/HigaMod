function J=get_J(nrb,u,n)
%% Ju for gauss-N
 nrb_der=nrbderiv(nrb);
 [pnt, jac] = nrbdeval (nrb, nrb_der, u);
Ju=zeros(3*numel(u),numel(u));
for i=1:numel(u)
    Ju(i*3-2:i*3,i)=jac(:,i);
end
%% Jp for gauss_n
Jp=zeros(3*numel(u),3*n);
for i=1:numel(u)
    for j=1:nrb.number;
    der(i,j)=nrbeval_der_p(nrb, j,u(i));
    Jp(i*3-2:i*3,j*3-2:j*3)=diag([der(i,j),der(i,j),der(i,j)]);
    end
 end
%% Jw fo G-N
Jw=zeros(3*numel(u),n);
for i=1:numel(u)
    for j=1:nrb.number;
    Jw(i*3-2:i*3,j)=nrbeval_der_w(nrb, j,u(i)) ;
    end
end
%% create J
J=[Ju Jp Jw];
end