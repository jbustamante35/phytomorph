function [bool] = isOnDisk(fileName,moniker)
    % init default
    %%%%%%%%%%%%%%%%%%%%%%%%%
    bool = 1;
    
    % if nargin == 1 then return
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 1;return;end
    
    % make sure that the fileName/path is long enough to search
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(fileName) > numel(moniker)
        bool = strcmp(fileName(1:numel(moniker)),moniker);
    else
        bool = 0;
    end
    
end