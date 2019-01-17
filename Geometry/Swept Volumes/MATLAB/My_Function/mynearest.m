function [xn yn zn index d]=mynearest(xb,yb,zb,x,y,z);
% given xb,yb,zb find the nearest point on x,y,z
xn=x(1);yn=y(1);zn=z(1);
for i=2:numel(x)
    if sqrt((xb-x(i))^2+(yb-y(i))^2+(zb-z(i))^2)<sqrt((xb-xn)^2+(yb-yn)^2+(zb-zn)^2)
        xn=x(i);yn=y(i);zn=z(i);
        d=sqrt((xb-x(i))^2+(yb-y(i))^2+(zb-z(i))^2);
        index=i;
    end
end
end
