function num=nearest_face(s,f,v);
dist=inf;
n=numel(f)/3;
for i=1:n
    p=v(f(i,:),:);
    xgp=sum(p(:,1))/3;
    ygp=sum(p(:,2))/3;
    zgp=sum(p(:,3))/3;
    if norm([xgp ygp zgp]-s)<dist
        num=i;
        dist=norm([xgp ygp zgp]-s);
    end
end
end
