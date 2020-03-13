function [sFileList] = fastOrderintoSets(FileList)
    
    parfor e = 1:numel(FileList)
        [PTH{e},NM{e},EXT{e}] = fileparts(FileList{e});
    end
    [uPTH] = unique(PTH);
    
    
    sFileList = cell(1,numel(uPTH));
    parfor u = 1:numel(uPTH)
        
        try
            
            idx = strcmp(PTH,uPTH{u});
            sFileList{u} = FileList(idx);
            tmp = NM(idx);
            nName = cellfun(@(X)str2num(X),tmp);
            [~,sidx] = sort(nName);
            sFileList{u} = sFileList{u}(sidx);
        catch
            sFileList{u} = {''};
        end
        
    end
    
end