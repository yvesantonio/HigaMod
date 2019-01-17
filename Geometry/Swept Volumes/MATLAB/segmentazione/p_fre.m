function p_fre(St,step)
if St.cont==1
    p_fre(St.next1,step);
    p_fre(St.next2,step);
end
plot_frenet(St.data(:,1),St.data(:,2),St.data(:,3),step)