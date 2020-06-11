function [bool] = isIRODS(fileName,moniker)
    if nargin == 1;moniker = '/iplant';end
    [bool] = isOnDisk(fileName,moniker);
end