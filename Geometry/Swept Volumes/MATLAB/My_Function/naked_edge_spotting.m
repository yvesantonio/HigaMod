clear all
close all
clc
[v,f,normals,name] = stlRead('C:\Users\Leo\Desktop\model.stl');
%stlPlot(v,f,name)
%% check if is naked edge (only belong to 1 face)
list=naked_edge(f);
%% chain what can be chained
l=1;
while list~=-1
    [chained_edge{l}(:) list]=chain(list);
    l=l+1;
end
%% plot naked edge chained togheter
stlPlot(v,f,name);hold on
for i=1:numel(chained_edge)
plot3(v(chained_edge{i},1),v(chained_edge{i},2),v(chained_edge{i},3),'Color',[1 0 i/numel(chained_edge)],'LineWidth',3);hold on
end
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
%% plot the dressed edge
mystlPlot(v,f,0.3);hold on
plot3(v(end-1:end,1),v(end-1:end,2),v(end-1:end,3),'rd','MarkerSize',3);hold on
