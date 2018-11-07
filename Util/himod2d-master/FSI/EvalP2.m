function uquad=EvalP2(u,FEM,ne,nqnx)

uquad=zeros(ne*nqnx,1);
for i=1:ne
    for q=1:nqnx
        uquad((i-1)*nqnx + q) = uquad((i-1)*nqnx + q)      ...
            + ( u(1 +(i-1)*2)*FEM(q,1)     ...
            +   u(2 +(i-1)*2)*FEM(q,2)     ...
            +   u(3 +(i-1)*2)*FEM(q,3) );
    end
end