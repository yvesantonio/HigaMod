%% script per plottare gli errori
clear all
close all
%Loading results of the simulation.
load caso8

% Setting grafica
C=100*L2(1,1);

% Errore L2
figure(1);
v=[1 2 4 8 16 32];
L2o=L2;

loglog(v,L2o,'-d','linewidth',3);
hold on;
%grid on

%loglog(mout,C*mout.^(-1),'--r','linewidth',3);%'O(1)',
loglog([1 32],C*[1 32].^(-2),'k','linewidth',2);
%loglog(mout,C*mout.^(-3),'--b','linewidth',3);%,'O(3)'
% loglog(mout,C*mout.^(-4),'--c','linewidth',3);%,'O(4)'
ylim([1e-6,C])


legend(['h=',h{1}],['h=',h{2}],['h=',h{3}],['h=',h{4}],['h=',h{5}],'Order=2');

set(gca,'FontUnits','points','FontSize',18)
xlabel('m')
title('L^2 - error')

%% Errore H1
clear all
close all
%Loading results of the simulation.
load caso8

figure
D=10*H1(1,1);
v=[1 2 4 8 16 32];
H1o=H1;
loglog(v,H1o,'-d','linewidth',3);
hold on;
%grid on
loglog(v,D*v.^(-1),'-k','linewidth',2);
%loglog(mout(:,1),D*mout(:,1).^(-2),'--g');
%loglog(mout(:,1),D*mout(:,1).^(-3),'--b');
legend(['h=',h{1}],['h=',h{2}],['h=',h{3}],['h=',h{4}],['h=',h{5}],'order=1','location','SouthWest');
ylim([1e-4,D])
set(gca,...
    'FontUnits','points', ...
    'FontSize',18);
xlabel('m')
title('H^1 - Error')
