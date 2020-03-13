function [varBlock] = functionHasVarblock(functionName)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find function
    % and open file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    funcPath = which(functionName);
    fileID = fopen(funcPath,'r');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check for variable block
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    textLine = fgetl(fileID);
    varBlock = false;
    varBlockLine = '<varBlock>';
    while (ischar(textLine) && ~varBlock)
        varBlock = contains(textLine,varBlockLine);
        textLine = fgetl(fileID);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % close file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(fileID);
end