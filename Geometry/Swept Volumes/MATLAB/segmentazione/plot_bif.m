function plot_bif(St,v)
if St.cont==1
mystlPlot_nofig(v,St.bif,1)
plot_bif(St.next1,v)
plot_bif(St.next2,v)
end
end