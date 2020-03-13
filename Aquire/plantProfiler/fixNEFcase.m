function [fileName] = fixNEFcase(fileName)
    if ~exist(fileName)
        fileName = strrep(fileName,'.nef','.NEF');
    end
end