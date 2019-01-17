function plotsection(St)
for i=1:numel(St.sec)
% plot3(St.sec{i}(:,1),St.sec{i}(:,2),St.sec{i}(:,3),'rd')
plot3(St.sec{i}(:,1),St.sec{i}(:,2),St.sec{i}(:,3),'r','LineWidth',3)
end
if St.cont==1
    plotsection(St.next1)
end

if St.cont==1
    plotsection(St.next2)
end
