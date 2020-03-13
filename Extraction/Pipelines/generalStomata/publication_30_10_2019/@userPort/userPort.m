classdef userPort < outPort
    
    methods
        
        function [obj] = userPort(oPath)
            if nargin == 0;oPath = './output/';end
            if isempty(oPath);oPath = './output/';end
            obj@outPort(oPath);
            %obj.initRemotePath();
        end
        
        function [] = initRemotePath(obj)
            
        end
        
        
    end
    

end