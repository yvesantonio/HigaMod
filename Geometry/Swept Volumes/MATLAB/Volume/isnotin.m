function is=isnotin(x,y,z,xx,yy,zz);
is=1;
for i=1:numel(xx)
    if x==xx(i) & y==yy(i) & z==zz(i)
        is=0;
    end
end
end
