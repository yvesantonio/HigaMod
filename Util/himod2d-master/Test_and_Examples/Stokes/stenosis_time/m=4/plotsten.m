close all
clear all
% data assimilation
t=[1,1.4,1.5,1.6,2,2.3];
n=6; 
figure
xlength=0.4;
ylength=0.2;
dx=0.5;
dy=0.33;
for k=1:n    
    pari = 1-mod(k,2);
    alt  = ceil(k/2)-1;
    h=subplot(3,2,k);
    set(h,'position',[0.05+pari*dx, 0.73-dy*alt, xlength, ylength]);
    importfile(['U_profile_x_fem_t=',num2str(t(k)),'.dat'])
    plot(data(:,2),data(:,1),'-','color',[139 136 120]/sum([139 136 120]),'linewidth',3);
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