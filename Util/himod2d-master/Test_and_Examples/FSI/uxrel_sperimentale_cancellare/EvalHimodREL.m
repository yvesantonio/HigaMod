function U=EvalHimodREL(Xhat,u,MB,FEM,fedof,nqnx)
ne=size(Xhat,2)/nqnx;
U=zeros(ne*nqnx,size(MB,1));
for k=3:size(MB,2)
    for j=1:size(MB,1)
        for i=1:ne
            for q=1:nqnx
                U((i-1)*nqnx + q,j) = U((i-1)*nqnx + q,j)        ...
                    + ( u( (k-1-2)*fedof +1+2 +(i-1)*2)*FEM(q,1)     ...
                    +   u( (k-1-2)*fedof +2+2 +(i-1)*2)*FEM(q,2)     ...
                    +   u( (k-1-2)*fedof +3+2 +(i-1)*2)*FEM(q,3) ) * ...
                    MB(j,k);
            end
        end
    end
end
for k=1:2
    for j=1:size(MB,1)
        i=1;
            for q=1:nqnx
                U((i-1)*nqnx + q,j) = U((i-1)*nqnx + q,j)        ...
                    + u(k +(i-1)*2)*FEM(q,1) * MB(j,k);
            end
    end
end