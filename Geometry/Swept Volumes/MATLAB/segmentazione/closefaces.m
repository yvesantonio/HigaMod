function [f v]=closefaces(f,v)
list=naked_edge(f);
l=1;
while list~=-1
    [chained_edge{l}(:) list]=chain(list);
    l=l+1;
end
% usare i baricentri per raccordare le centerline
%% create new faces from baricenter of linked edge 
for i=1:numel(chained_edge)
xb=sum(v(chained_edge{i},1))/numel(v(chained_edge{i},1));
yb=sum(v(chained_edge{i},2))/numel(v(chained_edge{i},2));
zb=sum(v(chained_edge{i},3))/numel(v(chained_edge{i},3));
v(end+1,:)=[xb yb zb];
for j=1:numel(chained_edge{i})-1
    f(end+1,:)=[chained_edge{i}(j:j+1) numel(v)/3];    
end
end
end