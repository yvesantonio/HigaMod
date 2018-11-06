function [] = MyWriteVTK(ffVTK,higaSOL)

    fileID = fopen(ffVTK);
    dataStruct = textscan(fileID,'%s %s %s %s %s');
    fclose(fileID);
    
    % POSITION OF THE INFORMATION ON THE NODES
    
    numbNodesPos = 5;
    aux = dataStruct{2}(numbNodesPos);
    auxx = aux{1};
    numbNodes = str2double(auxx);
    
    % POSITION OF THE INFORMATION ON THE CELLS
    
    numbCellsPos = numbNodesPos + numbNodes + 1;
    aux = dataStruct{2}(numbCellsPos);
    auxx = aux{1};
    numbCells = str2double(auxx);
    
    % POSITION OF THE INFORMATION ON THE CELL TYPES
    
    CellsTypePos = numbCellsPos + numbCells + 4;
    LookUpTabPos = CellsTypePos + 4;
    PntDataPos = LookUpTabPos + numbCells + 10;
    
    fileID = fopen(ffVTK);
    dataStruct2 = textscan(fileID,'%s','Delimiter','\n');
    fclose(fileID);
    
    % WRITE THE HIGAMOD SOLUTION .VTK FILE
    
    disp('STARTED WRITTING .VTK FILE')
    
    aux = dataStruct2{1};
    filename = 'HigaMod.vtk';
    fid = fopen(filename, 'w');
    for ii = 1:size(dataStruct2{1},1)
        auxx = aux{ii};
        
        if (ii < PntDataPos) || (ii >= PntDataPos + length(higaSOL))
            fprintf(fid, [auxx,'\n']);
        else
            fprintf(fid, [num2str(higaSOL(ii - PntDataPos + 1)),'\n']);
        end
        
    end
    fclose(fid);
    
    disp('FINISHED WRITTING .VTK FILE')

end