a = 1;
b = [1/4 1/4 1/4 1/4];
yx = filter(b,a,x);
yy=filter(b,a,y);
yz=filter(b,a,z);
yx=yx(numel(b):end);
yy=yy(numel(b):end);
yz=yz(numel(b):end);
figure
plot3(x,y,z,'--',yx,yy,yz,'-');grid on
legend('Original Data','Filtered Data')