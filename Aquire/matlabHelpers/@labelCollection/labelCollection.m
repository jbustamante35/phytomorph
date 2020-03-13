classdef labelCollection < doid
    properties
        name;
        description;
        labelSet;
    end
    
    methods
        function [obj] = labelCollection(name,description)
            obj.name = name;
            obj.description = description;
        end
        
        % 
        function [] = addLabel(obj,label)
            labelPtr = dptr.transformInput(label);
            obj.labelSet = [obj.labelSet labelPtr];
        end
    end
end