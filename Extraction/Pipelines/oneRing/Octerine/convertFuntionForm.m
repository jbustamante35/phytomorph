function [] = convertFuntionForm(functionName)




    funcPath = which(functionName);
    
    
    %[fList,pList] = matlab.codetools.requiredFilesAndProducts(funcPath);
    
    
    
    [fileID,errmsg] = fopen(funcPath,'r');
    textLine = fgetl(fileID);
    
    
    
    
    if contains(textLine,'function')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clip out the input variable names
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        var_start = strfind(textLine,'(');
        var_end = strfind(textLine,')');
        vars = textLine(var_start:var_end);
        vars = strrep(vars,'(',',');
        vars = strrep(vars,')',',');
        cidx = strfind(vars,',');
        for e = 1:(numel(cidx)-1)
            varName{e} = vars((cidx(e)+1):(cidx(e+1)-1));
        end
        
        
        
        varBlock = functionHasVarblock(functionName);
        writeVarblock(functionName);
        addVarToBlock(functionName,'test','This is the text for test.');
        
        
        
        if ~varBlock
            for e = 1:numel(varName)
                des{e} = num2str(e);
                %des{e} = inputdlg('Enter description for variable',varName{e},[1, dlgWidth]);
            end
        else
            
        end
        
        
        
        
        frewind(fileID)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over lines to see variale usage
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        varUse = [];
        textLine = fgetl(fileID);
        while ischar(textLine)
            boolVec = zeros(1,numel(varName),'logical');
            for e = 1:numel(varName)
                boolVec(e) = contains(textLine,varName{e});
            end
            varUse = [varUse;boolVec];
            textLine = fgetl(fileID);
        end
        
        
        
    end
    
end