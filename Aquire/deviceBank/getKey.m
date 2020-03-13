function [k] = getKey(kvp)
    k = '';
    try
        fidx = strfind(kvp,'_');
        k = kvp(2:(fidx(1)-1));
    catch
    end
end