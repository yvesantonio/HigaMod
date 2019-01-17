function [St,f,v]=add_mesh2(St,f,v)

    
if St.next1.cont==1
   [St.next1,f,v]=add_mesh2(St.next1,f,v);
end
if St.next2.cont==1
   [St.next2,f,v]=add_mesh2(St.next2,f,v);
end

St.lim=[];
for j=1:numel(St.sec)
    
[f v]=mesh_adapt(f,v,St.sec{j});

[x y z a]=nearest(St.sec{j}(1,1),St.sec{j}(1,2),St.sec{j}(1,3),v(:,1),v(:,2),v(:,3));
p=0;
p(1)=a;
for i=2:numel(St.sec{j})/3
    [x y z a]=nearest(St.sec{j}(i,1),St.sec{j}(i,2),St.sec{j}(i,3),v(:,1),v(:,2),v(:,3));
    if a~=p(end)
        p(end+1)=a;
    end
end
p=p';
for i=1:numel(p)-1
    limit(i,:)=[p(i) p(i+1)];
end
St.lim{end+1}=limit;
clearvars limit
end

disp('addato')






