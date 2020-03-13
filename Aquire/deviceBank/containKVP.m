function [valid] = containKVP(inString,kvpSet)
    valid = validPSON(inString);
    if ~valid;return;end
    for e = 1:numel(kvpSet)
        curKVP = kvpSet{e};
        matchKey = getKey(curKVP);
        matchValue = getValue(curKVP);
        valid = containsKeys(inString,{matchKey});
        if ~valid;break;end
        inKVP = findKVP(inString,matchKey);
        value = getValue(inKVP);
        if ~strcmp(matchValue,'*')
            valid = valid & strcmp(value,matchValue);
        end
        if ~valid;break;end
    end
end