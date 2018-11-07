set(gca,'FontSize',20)
plot(mediaestremi_plot(:,1),mediaestremi_plot(:,2),'linewidth',3)
hold on
plot([3,3],[1.65,2],'--k','linewidth',2)
plot([6,6],[1.65,2],'--k','linewidth',2)
ylim([1.65,2])