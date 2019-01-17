function [] = exportTxt(matrix,numbSamples,filterTime)

    name = inputname(1);
    fileName = [name,'_s',num2str(numbSamples),'_t',num2str(filterTime),'.txt'];

    fileID = fopen(fileName,'w');
    fprintf(fileID,'%.15f  \n',matrix);
    fprintf(fileID,'\n');
    fprintf(fileID,'%.15f  \n',numbSamples);
    fprintf(fileID,' \n');
    fprintf(fileID,'%.15f  \n',filterTime);
    fclose(fileID);

end