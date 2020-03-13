function [] = objectForm(object)
    fu.form = 'json'
    fu.isPtr = false;

    % file pointer
    if isa(object,'file')
        fu.form = 'file';
        fu.isPtr = true;
        return
    end
    % if json - could be pts
    if isa(object,'char')
        fu.form = 'json';
        tmp = jsondecode(object);
        fu.isPtr = isPtr(tmp);
    end
    
    
    
    if isa(object,'struct')
        fu.form = 'struct';
        fu.isPtr = isPtr(tmp);
        return
    end
end