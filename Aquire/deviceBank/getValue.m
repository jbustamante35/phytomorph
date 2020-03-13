function [v] = getValue(kvp)
    v = '';
    try
        fidx = strfind(kvp,'_');
        v = kvp((fidx(1)+1):end-1);
    catch
    end
end