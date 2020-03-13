function [kvp] = findKVP(inString,key)
    try
        kvp = '';
        if ~isempty(inString) && ~isempty(key)
            for e = 1:numel(inString)
                str(e) = -strcmp(inString(e),'{');
                stp(e) = strcmp(inString(e),'}');
            end
            S = cumsum(str+stp);
            stpIdx = find(S == 0);
            strIdx = find(str==-1 & S == -1);
            for e = 1:numel(strIdx)
                kvpList{e} = inString(strIdx(e):stpIdx(e));
                k = getKey(kvpList{e});
                if strcmp(k,key)
                    kvp = kvpList{e};
                    break
                end
            end
        end
    catch ME
        ME
    end
    
    %{
    kvp = '';
    try
        keyLength=numel(key);
        kidx = strfind(str,[key '_']);

        endIdx = strfind(str,'}');
        endIdx = endIdx(endIdx > kidx);
        endIdx = endIdx(1);


        strIdx = strfind(str,'{');
        strIdx = strIdx(strIdx < kidx);
        strIdx = strIdx(end);

        kvp = str(strIdx:endIdx);
    catch
    end
    %}
    
end