function x_=bound_apply(x_,m,n)
u=x_(1:m);
for i=1:numel(u)
    if x_(i)<0
        x_(i)=0;
    else
        if x_(i)>1
            x_(i)=0.9999999;
        end
    end
end
for i=0:n-1
    if x_(end-i)<0.05
        x_(end-i)=0.05;
    else
        if x_(end-i)>20
            x_(end-i)=20;
        end
    end
end
end