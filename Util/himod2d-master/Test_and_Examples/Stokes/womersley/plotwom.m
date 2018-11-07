close all
% data assimilation
t=[1,1.4,1.5,1.6,2,2.3];
n=6;
v=1:16;
v=2*v;
figure
for k=1:6
    
    subplot(3,2,k)
    x=linspace(0,1);
    plot(exactsol(x,t(k)),x,'linewidth',3,'color',[139 136 120]/sum([139 136 120]))
    hold on
    importfile(['U_profile_x=0.5t=',num2str(t(k)),'.dat'])
    plot(data(v,2),data(v,1),'o','MarkerEdgeColor','none','MarkerFaceColor','red');
    importfile(['lU_profile_x=0.5t=',num2str(t(k)),'.dat'])
    plot(data(v,2),data(v,1),'o','MarkerEdgeColor','none','MarkerFaceColor','blue');
    clear data
    grid on
    xlim([-0.1,0.1])
    title(['\fontsize{20}t=',num2str(t(k))])
    set(gca,'FontSize',20)
   
end