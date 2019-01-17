function [x, iter]=mygaus(A,b,x,tol,max_it)
D=diag(diag(A));
L=tril(A,-1);
U=triu(A,1);
B=-(D+L)\U;
g=(D+L)\b;
for i=1:max_it & abs((A*x-b))>tol
    x=B*x+g;
end
iter=i;
end
