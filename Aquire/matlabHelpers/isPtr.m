function [b] = isPtr(obj)
    b = strcmp(obj.type(1),'>');
end