function [tmpPath] = makeTempLocation()
    tmpPath = [tempname filesep];
    mkdir(tmpPath);
end