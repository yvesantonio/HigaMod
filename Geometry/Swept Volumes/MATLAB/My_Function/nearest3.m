function res=nearest3(pi,vertices)
N=numel(vertices)/3;
res=zeros(3,4);


d=inf;
for i=1:numel(vertices)/3
    dn=norm(pi-vertices(i,:));
    if dn<d
        res(1,:)=[vertices(i,:) i];
        d=dn;
    end
end
d=inf;
for i=1:numel(vertices)/3
    dn=norm(pi-vertices(i,:));
    is1=ismember(vertices(i,:),res(1,:));
    if dn<d & is1==0
        res(2,:)=[vertices(i,:) i];
        d=dn;
    end
end
d=inf;
for i=1:numel(vertices)/3
    dn=norm(pi-vertices(i,:));
    is1=ismember(vertices(i,:),res(1,:));
    is2=ismember(vertices(i,:),res(2,:));
    if dn<d & is1==0 & is2==0
        res(3,:)=[vertices(i,:) i];
        d=dn;
    end
end
end
    