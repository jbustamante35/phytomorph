function [b] = isa(this,type)
    %b = isa(this,type);
    b = builtin('isa',this,type);
    %b = strcmp(class(this),type);
end