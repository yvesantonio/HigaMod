function plot_all(St,v)
if St.cont==0
    if St.tuboexist==1
    nrbplot(St.VOL,[10 10 100])
    mystlPlot_nofig(v,St.tubo,0.7)
    end
else
if St.tuboexist==1
    mystlPlot_nofig(v,St.tubo,0.7)
    nrbplot(St.VOL,[50 50 500])
end
plot_all(St.next1,v)
plot_all(St.next2,v)
end

end
