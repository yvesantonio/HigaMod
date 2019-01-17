

%% crea struttura di backup con le inlet ed outlet section chiuse
f=faces;
v=vertices;
[f v]=closefaces(f,v);
mystlPlot(v,f,1)
%% creo le sezioni interne da distaccare, uso ricorsione
[St f v]=create_section(St,v,f); %
f=faces;
v=vertices;

%% modifica le mesh
[St f v]=add_mesh(St,f,v);
%% disegna quello che si è ottenuto modificando le mesh
mystlPlot(v, f, 0.5);hold on
plotlimit(St,v)
plotc(St)
plotsection(St)
pause(1)
%% inizializza e separa i vri pezzi mettendo le facce che li costituiscono nella St
St.tuboexist=1;
St.start=1;
St=segmentazione(St,f,v);

%% plot segmented structure
figure;hold on
plot_segment(St,v)
%% 
