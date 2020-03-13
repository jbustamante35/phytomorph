classdef containerMsg < doid
    
    
    properties
        data;
    end
    
    
    methods
        function [obj] = containerMsg(data)
            obj.data = data;
        end
    end
    
end