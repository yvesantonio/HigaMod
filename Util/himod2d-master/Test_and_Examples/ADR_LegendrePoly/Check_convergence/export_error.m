load caso10

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
