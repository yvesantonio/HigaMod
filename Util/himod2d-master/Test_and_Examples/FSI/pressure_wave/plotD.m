clc
clear all
close all

v=0.4:0.4:1.6;
C=hsv(length(v));
for i=1:length(v)
    filename=['DU=',num2str(v(i)),'.dat'];
    importfile(filename);
    plot(data(:,1),data(:,2),'color',C(i,:),'linewidth',3)
    hold on
    clear data
    filename=['DD=',num2str(v(i)),'.dat'];
    importfile(filename);
    plot(data(:,1),data(:,2),'color',C(i,:),'linewidth',3)
    clear data
end