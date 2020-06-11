function [localList] = xfer_get_p(remoteList,localPath)

    mmkdir(localPath);

    
    parfor e = 1:numel(remoteList)
        localList{e} = xfer_get(remoteList{e},localPath);
    end
    
end


