function [P,Pd]=polynomial_legendre(m)
P=cell(m,1);
Pd=cell(m,1);
P{1} = 1;
Pd{1}=polyder(P{1});
P{2} = [1 0 ];
Pd{2}=polyder(P{2});
for n=1:m-2
    P{n+2}=1/(n+1)*((2*n+1)*[P{n+1} 0]-n*[0 0 P{n}]);
    Pd{n+2}=polyder(P{n+2});
end

for n=1:m
    C=sqrt((2*(n-1)+1)/2);
    P{n}=P{n}*C;
    Pd{n}=Pd{n}*C;
end
