function N=LS_N_point(x,y,z,d,n)

u=linspace(0,1,numel(x)); u(end)=0.99999;
m=numel(u);
U=zeros(1,n+d+1); % automatic node creation
U(end-d:end)=1;
mid=numel(U)-2*d-2;
for i=1:mid
    U(i+d+1)=1/(mid+1)*i;
end
%%
for i=1:n
    for j=1:numel(u)
    element(j,i)=Nip(i,u(j),d,U);
    end
end
%%
bx=x';
by=y';
bz=z';
xL=element\bx;
xL(1,1)=x(1);
xL(end,1)=x(end);
yL=element\by;
yL(1,1)=y(1);
yL(end,1)=y(end);
zL=element\bz;
zL(1,1)=z(1);
zL(end,1)=z(end);
pL(1:3,1)=[x(1);y(1);z(1)];
for i=1:numel(xL)
    pL(i*3-2:i*3,1)=[xL(i);yL(i);zL(i)];
end
%% nurbs toolbox prova
N=nrbmak([xL';yL';zL'],U);

end