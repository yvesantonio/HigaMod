clear all
close all
DELIMITER = ' ';
HEADERLINES = 1;
newData1 = importdata('U_profile_x=2t=100.dat', DELIMITER, HEADERLINES);
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
plot(5-5*data(:,1).^2,data(:,1),'color','red','linewidth',3)
hold on
plot(data(:,2),data(:,1),'d','MarkerFaceColor','blue');
xlim([-1,6])
grid on
dim=20;
% legend('Exact profile','HiMod solution','Location','SouthEast')
% title('$U_x$ velocity, x=2','interpreter','latex','FontSize',dim)
ylabel('y','FontSize',dim)
% xlabel('$U_x$','interpreter','latex','FontSize',dim)
set(gca,'FontSize',dim)

