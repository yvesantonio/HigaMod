function [B1L,B2L,B1R,B2R] = createMatrices(b1,b2,n1,n2)
B1L=zeros(n2,n1);
B2L=zeros(n2,n1);
B1R=zeros(n2,n1);
B2R=zeros(n2,n1);
N=n1*n2;
for k=1:n1 %y
    for j=1:n2 %x
        B1R(j,k)=b1(j+(k-1)*n2);
        B2R(j,k)=b2(j+(k-1)*n2);
    end
end
for k=1:n1
    for j=1:n2
        B1L(j,k)=b1(N+j+(k-1)*n2);
        B2L(j,k)=b2(N+j+(k-1)*n2);
    end
end
end

