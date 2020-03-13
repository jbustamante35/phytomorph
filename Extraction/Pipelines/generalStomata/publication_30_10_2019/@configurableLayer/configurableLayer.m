classdef (Abstract) configurableLayer < generalizedFunctionLayer

    properties
        isConfigured;
    end

    methods
        function [obj] = configurableLayer(name)
            obj = obj@generalizedFunctionLayer(name);
            obj.isConfigured = false;
        end
        
       
    end
    
    methods (Abstract)
        configure(varargin);
    end
    
end