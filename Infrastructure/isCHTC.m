function [bool] = isCHTC(fileName,moniker)
    if nargin == 1;moniker = '/chtc';end
    [bool] = isOnDisk(fileName,moniker);
end