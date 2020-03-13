function [b] = comparePtrs(ptrA,ptrB)
    if isPtr(ptrA);ptrA = ptrA.refs.uuid;else;ptrA=ptrA.uuid;end
    if isPtr(ptrB);ptrB = ptrB.refs.uuid;else;ptrB=ptrB.uuid;end
    b = strcmp(ptrA,ptrB);
end