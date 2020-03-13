classdef stringBuffer < handle

    properties (SetObservable = true)
        string;
    end
    
    
    methods
        
        function [obj] = stringBuffer()
            obj.string = '';
        end
        
        % hard coded to call processCommand
        function [] = attachListener(obj,listener)
            addlistener(obj,'string','PostSet',@(src,evnt)listener.processCommand(src,evnt)); % Add ob
        end
    end
end