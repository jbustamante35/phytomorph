classdef resultCollection < dataCollection
    
    properties
        fsStatus;
    end
    
    methods
        
        function [obj] = resultCollection(fileList)
            obj@dataCollection(fileList);
        end
        
        function [] = setStatus(this,status)
            this.fsStatus = status;
        end
        
    end
    
end