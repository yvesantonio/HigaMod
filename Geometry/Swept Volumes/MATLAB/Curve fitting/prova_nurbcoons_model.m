clc
clear all
close all
addpath('C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\NURBS',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\nurbs_toolbox',...
'C:\Users\Leo\Documents\MATLAB\Add-Ons\Collections\My_Function');
%% LS con curva di grado 3 e 6 punti
load('data_model.mat')
x=data_model(:,1); y=data_model(:,2); z=data_model(:,3);
figure; plot3(x,y,z,'r'); grid on; title('curva iniziale r' )
[pti  Cpti]=nrb_approx(x,y,z,11,2,1.e-7,20,100);
hold on
%plot3(pti(1,:),pti(2,:),pti(3,:),'b')
%  hold on;
% figure
% nrbplot(Cpti,1000)
%%
sud=round(numel(x)/4);
x1=x(1:sud+1)';
y1=y(1:sud+1)';
z1=z(1:sud+1)';
x4=x(sud+1:2*sud)';
y4=y(sud+1:2*sud)';
z4=z(sud+1:2*sud)';
x3=x(end:-1:end-sud)';
y3=y(end:-1:end-sud)';
z3=z(end:-1:end-sud)';
x2=x(end-sud:-1:2*sud)';
y2=y(end-sud:-1:2*sud)';
z2=z(end-sud:-1:2*sud)';
figure
n=10;
d=2;

[N1pti  U1]=nrb_approx(x1',y1',z1',n,d,1.e-7,20,100);
hold on;
plot3(N1pti(1,:),N1pti(2,:),N1pti(3,:),'r','LineWidth',2)

[N2pti  U2]=nrb_approx(x2',y2',z2',n,d,1.e-7,20,100);
hold on;
plot3(N2pti(1,:),N2pti(2,:),N2pti(3,:),'g','LineWidth',2)

[N3pti  V1]=nrb_approx(x3',y3',z3',n,d,1.e-7,20,100);
hold on;
plot3(N3pti(1,:),N3pti(2,:),N3pti(3,:),'b','LineWidth',2)

[N4pti  V2]=nrb_approx(x4',y4',z4',n,d,1.e-7,20,100);
hold on;
plot3(N4pti(1,:),N4pti(2,:),N4pti(3,:),'y','LineWidth',2)

%figure
nrbeval (U1, U1.knots(1)) - nrbeval (V1, V1.knots(1))
nrbeval (U1, U1.knots(end)) - nrbeval (V2, V2.knots(1))
nrbeval (U2, U2.knots(1)) - nrbeval (V1, V1.knots(end))
nrbeval (U2, U2.knots(end)) - nrbeval (V2, V2.knots(end))
P = nrbcoons(U1, U2, V1, V2);
nrbplot(P,[20 20]);hold on;
plot3(x,y,z,'b')