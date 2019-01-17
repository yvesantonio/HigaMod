function plot_vol_and_segment(St,v)
if St.cont==1
    plot_vol_and_segment(St.next1,v);
    plot_vol_and_segment(St.next2,v);
end  
if St.tuboexist==1
figure
nrbplot(St.VOL,[50 50 500])    
mystlPlot_nofig(v,St.tubo,0.75)
end
end