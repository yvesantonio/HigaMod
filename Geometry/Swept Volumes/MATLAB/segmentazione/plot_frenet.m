function plot_frenet(x,y,z,step);

[Tg Nm Bn]=myfrenet(x,y,z);
for i=1:numel(x)-1
    if Tg(i,:)*Tg(i+1,:)'<0
        Tg(i+1,:)=-Tg(i+1,:);
    end
end
%% plot vector tg nm bn

% % DEBUG: PLOT THE VECTOR TANGENT TO THE CENTERLINE IN FIVE DIFFERENT
% % LOCATIONS
hold on;
color=[randi([0,10])/10 randi([0,10])/10 randi([0,10])/10];
for i = 1:step:numel(x)
    pInit=[x(i); y(i); z(i)];
    pFinal=[x(i)+Tg(i,1); y(i)+Tg(i,2); z(i)+Tg(i,3)]; hold on     
    vectarrow(pInit,pFinal)   %plot3([pInit(1) pFinal(1)],[pInit(2) pFinal(2)],[pInit(3) pFinal(3)],'LineWidth',1,'Color',color)
end
color=[randi([0,10])/10 randi([0,10])/10 randi([0,10])/10];
for i = 1:step:numel(x)
    pInit=[x(i); y(i); z(i)];
    pFinal=[x(i)+Nm(i,1); y(i)+Nm(i,2); z(i)+Nm(i,3)];hold on      
    vectarrow(pInit,pFinal) %plot3([pInit(1) pFinal(1)],[pInit(2) pFinal(2)],[pInit(3) pFinal(3)],'LineWidth',1,'Color',color)
end
color=[randi([0,10])/10 randi([0,10])/10 randi([0,10])/10];
for i = 1:step:numel(x)
    pInit=[x(i); y(i); z(i)];
    pFinal=[x(i)+Bn(i,1); y(i)+Bn(i,2); z(i)+Bn(i,3)];  hold on    
    vectarrow(pInit,pFinal) %plot3([pInit(1) pFinal(1)],[pInit(2) pFinal(2)],[pInit(3) pFinal(3)],'LineWidth',1,'Color',color)
end
