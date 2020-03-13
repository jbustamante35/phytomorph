function [] = writeVarblock(functionName)
    has_varBlock = functionHasVarblock(functionName);
    if ~has_varBlock
        funcPath = which(functionName);
        fileID = fopen(funcPath,'at');

        commentBlockStart = '%{';
        commentBlockStop = '%}';


        varBlockStart = '<varBlock>';
        varBlockStop = '</varBlock>';

        fprintf(fileID,'\n\n');


        fprintf(fileID,'\n');
        fprintf(fileID,'%s',commentBlockStart);
        fprintf(fileID,'\n');
        fprintf(fileID,'%s',varBlockStart);
        fprintf(fileID,'\n');


        fprintf(fileID,'\n');
        fprintf(fileID,'%s',varBlockStop);
        fprintf(fileID,'\n');
        fprintf(fileID,'%s',commentBlockStop);
        fprintf(fileID,'\n');


        fclose(fileID);
    end
end