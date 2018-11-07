clear all
close all
DELIMITER = ' ';
HEADERLINES = 1;
newData1 = importdata('Pressure_t=100.dat', DELIMITER, HEADERLINES);
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
plot(data(:,1),5-data(:,1),'color','red','linewidth',3)
hold on
plot(data(:,1),data(:,2),'d','MarkerFaceColor','blue');
xlim([0,5])
grid on
dim=20;
% legend('Exact profile','HiMod solution','Location','SouthEast')
% title('$U_x$ velocity, x=2','interpreter','latex','FontSize',dim)
xlabel('x','FontSize',dim)
% xlabel('$U_x$','interpreter','latex','FontSize',dim)
set(gca,'FontSize',dim)

