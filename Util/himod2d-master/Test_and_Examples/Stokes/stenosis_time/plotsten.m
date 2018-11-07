close all
% data assimilation
t=[1,1.4,1.5,1.6,2,2.3];
n=6; 
figure
for k=1:n    
    subplot(3,2,k)
    importfile(['U_profile_x_fem_t=',num2str(t(k)),'.dat'])
    plot(data(:,2),data(:,1),'-k','linewidth',3);
    clear data
    hold on  
    importfile(['PFU_profile_x=4t=',num2str(t(k)),'.dat'])
    plot(data(:,2),data(:,1),'b','linewidth',3);
    clear data
    importfile(['ExU_profile_x=4t=',num2str(t(k)),'.dat'])
    plot(data(:,2),data(:,1),'r','linewidth',3);
    clear data
%     importfile(['NewU_profile_x=4t=',num2str(t(k)),'.dat'])
%     plot(data(:,2),data(:,1),'o','MarkerEdgeColor','none','MarkerFaceColor','red');
%     clear data
    grid on
    xlim([-0.0015,0.003])
    set(gca,'XTick',[-0.0015,0,0.001,0.0025])
    set(gca,'XTickLabel',{'-0.0015','0','0.001','0.0025'})
    title(['\fontsize{20}t=',num2str(t(k))])
    set(gca,'FontSize',20)
end