function [bool] = isSQUID(fileName,moniker)
    if nargin == 1;moniker = '/squid/ndmiller';end
    [bool] = isOnDisk(fileName,moniker);
end