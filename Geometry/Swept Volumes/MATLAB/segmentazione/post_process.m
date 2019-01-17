function St=post_process(St,v)
ne=nearest_naked_edge(St.tubo,v,St.data(St.start,:));
nf=[];
[nf f]=segmento(nf,St.tubo,ne(1:2),ne(3),St,v);
St.tubo=newf(f);
disp('finito_post_proc');
%mystlPlot(v, St.tubo, 0.3);hold on
