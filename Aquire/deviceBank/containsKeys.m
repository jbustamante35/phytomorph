function [ret] = containsKeys(str,keys)
    ret = true;
    for e = 1:numel(keys)
        tmp = findKVP(str,keys{e});
        ret = ~isempty(tmp)&ret;
    end
end