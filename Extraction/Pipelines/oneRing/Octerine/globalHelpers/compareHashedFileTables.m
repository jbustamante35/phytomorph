function [CList,Cindex] = compareHashedFileTables(A,B)

    AKeys = A.hashValue;
    BKeys = B.hashValue;

    Ckeys = setdiff(AKeys,BKeys,'stable');


    for e = 1:numel(Ckeys)
        Cindex(e) = find(strcmp(AKeys,Ckeys{e}));
        CList{e} = A.wholeName(Cindex(e));
    end

end