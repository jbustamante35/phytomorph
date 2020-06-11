function [FileList] = ldig(FilePath,FileList,FileExt,verbose)
    try
        
        
        
        %locate -d /mnt/scratch4/spaldingdata.db --regex "^/mnt/spaldingdata/Drone_Imagery/Arlington_2018/.*(\.JPG|\.jpg)$" | grep jpg
        dbList{1} = '/mnt/scratch4/spaldingdata.db';
        dbList{2} = '/mnt/scratch4/tetra.db';
        dbList{3} = '/mnt/scratch4/spaldingimages.db';
        dbList{4} = '/mnt/scratch4/snapper.db';
        
        
        % build up the command string
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CMD = ['locate '];
        for d = 1:numel(dbList)
            CMD = [CMD ' -d ' dbList{d} ' '];
        end
        CMD = [CMD '--regex '];
        CMD = [CMD '"^' FilePath '.*'];

        
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
        end
        
        

        [r,o] = system(CMD);
        oidx = strfind(o,char(10));
        oidx = [0 oidx];
        for e = 1:(numel(oidx)-1)
            FileList{e} = o(oidx(e)+1:oidx(e+1)-1);
        end

        
    catch ME
        ME
    end
    
end

