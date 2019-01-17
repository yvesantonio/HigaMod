function plotc(C)

if C.cont==0
    plot3(C.data(:,1),C.data(:,2),C.data(:,3),'Color',[randi([0,10])/10 randi([0,10])/10 randi([0,10])/10],'LineWidth',2)
else
    plotc(C.next1)
    plotc(C.next2)
    plot3(C.data(:,1),C.data(:,2),C.data(:,3),'Color',[randi([0,10])/10 randi([0,10])/10 randi([0,10])/10],'LineWidth',2)
end
    