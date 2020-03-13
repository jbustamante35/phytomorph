function [] = addVarToBlock(functionName,varName,varDes)
    has_varBlock = functionHasVarblock(functionName);
    
    if has_varBlock
        funcPath = which(functionName);
        [pth,nm,ext] = fileparts(funcPath);
        newFile = [pth filesep nm '_tmp' ext];
        fileID = fopen(funcPath,'r');
        fileID_new = fopen(newFile,'w');
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read and copy upto the varBlock
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        textLine = fgets(fileID);
        varBlock = false;
        varBlockLine = '<varBlock>';
        fileStore_PRE = {};
        while (ischar(textLine) && ~varBlock)
            % store the text data
            fileStore_PRE{end+1} = textLine;
            varBlock = contains(textLine,varBlockLine);
            % copy the text data into the new file
            fprintf(fileID_new,'%s',textLine);
            % get next line
            textLine = fgets(fileID);
        end
        % copy the text data into the new file
        fprintf(fileID_new,'%s',textLine);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read and copy upto the varBlock
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
        varStart = '<var name=''#NAME#'' description=''#DES#''>';
        varStop = '</var>';
        
        %fprintf(fileID,'\n');
        varLine = strrep(varStart,'#NAME#',varName);
        varLine = strrep(varLine,'#DES#',varDes);
        
        fprintf(fileID_new,'%s',varLine);
        fprintf(fileID_new,'\n');
        fprintf(fileID_new,varName);
        fprintf(fileID_new,'\n');
        fprintf(fileID_new,'%s',varStop);
        fprintf(fileID_new,'\n');
        
        
        
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read and store after the required location
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fileStore_POST = {};
        while (ischar(textLine))
            fileStore_POST{end+1} = textLine;
            fprintf(fileID_new,'%s',textLine);
            textLine = fgets(fileID);
        end
        
        

        fclose(fileID);
        fclose(fileID_new);
    end
end