function [bool] = isLOCAL(fileName)
    bool = ~isCHTC(fileName) && ~isIRODS(fileName) && ~isCOLD(fileName) && ~isSQUID(fileName);
end