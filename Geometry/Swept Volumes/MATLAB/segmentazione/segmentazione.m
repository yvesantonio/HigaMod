function St=segmentazione(St,f,v)
St.tubo=[];
if St.tuboexist
ne=nearest_naked_edge(f,v,St.data(St.start,:));
nf=[];
[nf f]=segmento(nf,f,ne(1:2),ne(3),St,v);
St.tubo=nf;
f=newf(f);
disp('finito');
mystlPlot(v, St.tubo, 0.3);hold on
pause(1)
end
if St.cont==1
ne=nearest_naked_edge(f,v,St.data(end,:));
nf=[];
[nf f]=segmento(nf,f,ne(1:2),ne(3),St,v);
St.bif=nf;
f=newf(f);
disp('finito')
mystlPlot(v, St.bif, 0.3);hold on
pause(1)
St.next1=segmentazione(St.next1,f,v);
St.next2=segmentazione(St.next2,f,v);
 
end
