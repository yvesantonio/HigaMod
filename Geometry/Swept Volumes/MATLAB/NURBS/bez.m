function C=bez(P)
u=0:0.01:1;
n=numel(P(:,1))-1;
    for i=0:n
        B(i+1,:)=factorial(n)/factorial(i)/factorial(n-i)*u.^i.*(1-u).^(n-i);
    end
    C=zeros(3,numel(u));
% Y=f(X)
for r=1:3 % x y z
for k=1:n+1 % tutti i punti
    C(r,:)=C(r,:)+P(k,r)*B(k,:);
end
end
end