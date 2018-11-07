close all
% data assimilation
t=[1,1.4,1.5,1.6,2,2.3];
n=6;
v=1:16;
v=2*v;
figure
for k=1:6
    importfile(['U_profile_x=0.5t=',num2str(t(k)),'.dat'])
    subplot(3,2,k)
    x=linspace(0,1);
    plot(exactsol(x,t(k)),x,'linewidth',3,'color','red')
    hold on
    plot(data(v,2),data(v,1),'d','MarkerFaceColor','blue');
    grid on
    xlim([-0.1,0.1])
    title(['\fontsize{20}t=',num2str(t(k))])
    set(gca,'FontSize',20)
    clear data
end