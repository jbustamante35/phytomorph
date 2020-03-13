function [inStruct] = addToStruct(inString,inStruct)
    try
        if ~isempty(inString)
            for e = 1:numel(inString)
                str(e) = -strcmp(inString(e),'{');
                stp(e) = strcmp(inString(e),'}');
            end
            S = cumsum(str+stp);
            stpIdx = find(S == 0);
            strIdx = find(str==-1 & S == -1);
            for e = 1:numel(strIdx)
                kvpList{e} = inString(strIdx(e):stpIdx(e));
                keys{e} = getKey(kvpList{e});
                values{e} = getValue(kvpList{e});
                if isfield(inStruct,keys{e})
                    inStruct.(keys{e}) = [inStruct.(keys{e}) values{e}];
                else
                    inStruct.(keys{e}) = values{e};
                end
            end
        end
    catch ME
        ME
    end
end