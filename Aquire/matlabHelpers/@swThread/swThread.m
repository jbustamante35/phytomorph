classdef swThread < wThread
    % cluster of worker threads
    properties
        
        
    end
    
    
    methods
        
        function [obj] = swThread(sharedToMain,clusterN)
            % super constructor
            obj@wThread(sharedToMain);
        end
        
        function [] = attachOutChannel(obj,entry)
            here = 1;
        end
         
    end
    
end