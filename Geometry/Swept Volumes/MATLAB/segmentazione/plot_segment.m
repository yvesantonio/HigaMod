function plot_all(St,v)
if St.cont==0 & St.tubo_exist==1
    nrbplot(St.tubo,[20 20 200])
else
mystlPlot_nofig(v,St.bif,1)
end
