load caso4
% Exporter per tikz
for k=1:NM
    L2(k,NH+1)=L2(1,1)/(2^k);
    L2(k,NH+2)=L2(1,1)/(4^k);
    L2(k,NH+3)=L2(1,1)/(8^k);
    L2(k,NH+4)=L2(1,1)/(16^k);
    H1(k,NH+1)=H1(1,1)/(2^k);
    H1(k,NH+2)=H1(1,1)/(4^k);
    H1(k,NH+3)=H1(1,1)/(8^k);
    H1(k,NH+4)=H1(1,1)/(16^k);
end
fil=fopen(['Caso',num2str(test),'_L2.dat'],'w');
fprintf(fil,'#h= ');
for i=1:NH-1
    fprintf(fil,[h{i},',']);
end
fprintf(fil,[h{NH},'\n']);
for k=1:NM
    fprintf(fil,'%d, ',mout(k,1));
    for i=1:NH+3
        fprintf(fil,'%12.16f, ',L2(k,i));
    end
    fprintf(fil,'%12.16f\n',L2(k,NH+4));
end
fprintf(fil,'#-----\n');
fclose(fil);

fil=fopen(['Caso',num2str(test),'_H1.dat'],'w');
fprintf(fil,'#h= ');
for i=1:NH-1
    fprintf(fil,[h{i},',']);
end
fprintf(fil,[h{NH},'\n']);
for k=1:NM
    fprintf(fil,'%d, ',mout(k,1));
    for i=1:NH+3
        fprintf(fil,'%12.16f, ',H1(k,i));
    end
    fprintf(fil,'%12.16f\n',H1(k,NH+4));
end
fprintf(fil,'#-----\n');
fclose(fil);


dim=size(L2);

righe=6;
colonne=5;
fprintf('\\begin{table}\n')
fprintf('\\centering\n')
fprintf('\\begin{tabular}{cccccc}\n')

fprintf('\\toprule\n')
fprintf('m&h=0.1&h=0.05&h=0.025&h=0.0125&h=0.00625\\\\\n')
fprintf('\\midrule\n')
for i=1:righe
	fprintf('%d  ', 2^(i-1));
	for j=1:colonne
		fprintf('& %1.2d  ', L2(i,j));
	end
	fprintf('\\\\ \n');
end
fprintf('\\bottomrule\n')
fprintf('\\end{tabular}\n')
fprintf('\\caption{Case test ???? -- error in $L^2$}\n');
fprintf('\\end{table}\n')


fprintf('\\begin{table}\n')
fprintf('\\centering\n')
fprintf('\\begin{tabular}{cccccc}\n')

fprintf('\\toprule\n')
fprintf('m&h=0.1&h=0.05&h=0.025&h=0.0125&h=0.00625\\\\\n')
fprintf('\\midrule\n')
for i=1:righe
	fprintf('%d  ', 2^(i-1));
	for j=1:colonne
		fprintf('& %.2d  ', H1(i,j));
	end
	fprintf('\\\\ \n');
end
fprintf('\\bottomrule\n')
fprintf('\\end{tabular}\n')
fprintf('\\caption{Case test ???? -- error in $H^1$}\n');
fprintf('\\end{table}\n')
