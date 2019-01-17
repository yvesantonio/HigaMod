function plotlimit(St,vertices)

if St.cont==1
    plotlimit(St.next1,vertices);
end

if St.cont==1
    plotlimit(St.next2,vertices);
end

for i=1:numel(St.lim)   
plot3(vertices(St.lim{i}(:,1),1),vertices(St.lim{i}(:,1),2),vertices(St.lim{i}(:,1),3),'LineWidth',3)
end
