function [bool] = isCOLD(fileName,moniker)
    if nargin == 1;moniker = '/mnt/spaldingarchive';end
    [bool] = isOnDisk(fileName,moniker);
end