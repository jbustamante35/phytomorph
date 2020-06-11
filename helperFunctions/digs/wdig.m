function [FileList] = wdig(FilePath,FileList,FileExt,verbose)
    try
        % build up the command string
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CMD = ['mc --json find ' FilePath(2:end) ' '];
        CMD = [CMD '--regex ".*'];
        
        if ~isempty(FileExt)
            % build up the extString
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            extString = '';
            for e = 1:numel(FileExt)
                extString = [extString '|\.' FileExt{e}]; 
            end
            extString(1) = '(';
            extString = [extString ')'];
            CMD = [CMD extString '$"'];
        else
            CMD = [CMD '"'];
        end
        
        
        [o,result] = system(CMD);
        result = strtrim(result);
        fidx1 = strfind(result,'{');
        fidx2 = strfind(result,'}');
        for e = 1:numel(fidx1)
            snip = result(fidx1(e):fidx2(e));
            fileStruct(e) = jsondecode(snip);
            FileList{end+1} = ['/' fileStruct(e).key];
        end
       
        
        
    catch ME
        getReport(ME)
    end
    
end
