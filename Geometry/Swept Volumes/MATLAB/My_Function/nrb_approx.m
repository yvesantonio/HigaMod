function [  approx  N  ]=nrb_approx(x,y,z,n,d,tol,max_iter_i,max_ieter_lsearch);
% x,y,z = sample of what to interpolate, must be column vector
% n = number of control points
% d = degree of spline
% num_pt_eval = number of eq. spaced point on which the best approx will be
% evaluated
% tol = tolerance for stop the while loop
% max_iter_k,max_ieter_k_l = limit for iteration
u=linspace(0,1,numel(x)); u(end)=0.99999999;
m=numel(u);
U=zeros(1,n+d+1); % automatic node creation
U(end-d:end)=1;
mid=numel(U)-2*d-2;
for i=1:mid
    U(i+d+1)=1/(mid+1)*i;
end
% create the system to find the first approximazion to put into the main
% script
clearvars element
for i=1:n
    for j=1:m
    element(j,i)=Nip(i,u(j),d,U);
    end
end
% x1,y1,z1 fist approx for control points
x1=element\x;
y1=element\y;
z1=element\z;
% set beginning and end at bound of curve
x1(1)=x(1); x1(n)=x(end);
y1(1)=y(1); y1(n)=y(end);
z1(1)=z(1); z1(n)=z(end);
[approx N]=g_n_opt(u,x1,y1,z1,x,y,z,U,tol,max_iter_i,max_ieter_lsearch);
end

