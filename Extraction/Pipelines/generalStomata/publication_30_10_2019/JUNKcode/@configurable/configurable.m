classdef (Abstract) configurableLayer

    properties
        isConfigured;
    end

    
    methods
        [c] = isConfigured()
    end
    
    methods (Abstract)
        configure(obj,varargin);
        
    end
    
end