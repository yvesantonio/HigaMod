function [pie t r]=point_in_edge(p,v1,v2)
pie=0;

t=abs((p-v1)./(v2-v1));

if t(1)<1 & abs(t(1)-t(2))<1.e-2 & abs(t(1)-t(3))<1.e-2
    pie=1;
end
end