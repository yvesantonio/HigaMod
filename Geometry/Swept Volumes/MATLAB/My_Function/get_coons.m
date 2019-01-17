function P=get_coons(surf,n,d)
%%
x=surf(:,1); y=surf(:,2); z=surf(:,3);
sud=round(numel(x)/4);
x1=x(1:sud+1)';
y1=y(1:sud+1)';
z1=z(1:sud+1)';
x4=x(sud+1:2*sud)';
y4=y(sud+1:2*sud)';
z4=z(sud+1:2*sud)';
x3=x(end:-1:end-sud)';
y3=y(end:-1:end-sud)';
z3=z(end:-1:end-sud)';
x2=x(end-sud:-1:2*sud)';
y2=y(end-sud:-1:2*sud)';
z2=z(end-sud:-1:2*sud)';
[N1pti  U1]=nrb_approx(x1',y1',z1',n,d,1.e-7,20,100);
[N2pti  U2]=nrb_approx(x2',y2',z2',n,d,1.e-7,20,100);
[N3pti  V1]=nrb_approx(x3',y3',z3',n,d,1.e-7,20,100);
[N4pti  V2]=nrb_approx(x4',y4',z4',n,d,1.e-7,20,100);
P = nrbcoons(U1, U2, V1, V2);
end