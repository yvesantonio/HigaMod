function U=EvalHimod(Xhat,u,MB,FEM,fedof,nqnx)
ne=size(Xhat,2)/nqnx;
U=zeros(ne*nqnx,size(MB,1));
for k=1:size(MB,2)
    for j=1:size(MB,1)
        for i=1:ne
            for q=1:nqnx
                U((i-1)*nqnx + q,j) = U((i-1)*nqnx + q,j)        ...
                    + ( u( (k-1)*fedof +1 +(i-1)*2)*FEM(q,1)     ...
                    +   u( (k-1)*fedof +2 +(i-1)*2)*FEM(q,2)     ...
                    +   u( (k-1)*fedof +3 +(i-1)*2)*FEM(q,3) ) * ...
                    MB(j,k);
            end
        end
    end
end

