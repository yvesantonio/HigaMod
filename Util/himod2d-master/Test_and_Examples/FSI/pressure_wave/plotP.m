clc
clear all
close all

v=0.2:0.4:2.4;
C=hsv(length(v));
for i=1:length(v)
    filename=['Pressure_t=',num2str(v(i)),'.dat'];
    importfile(filename);
    plot(data(:,1),data(:,2),'color',C(i,:),'linewidth',3)
    hold on
end