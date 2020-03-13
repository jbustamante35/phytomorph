classdef stateCube_lite < doid
    
    properties
    
         
        type;
        startFunction;
        channelGenFunc;
        pauseDuration;
    
    end
    
    
    methods
        
        function [obj] = stateCube_lite(type,startFunction,channelGenerationFunction,pauseDuration)
            obj.type = type;
            obj.startFunction = startFunction;
            obj.channelGenerationFunction = channelGenerationFunction;
            obj.pauseDuration = pauseDuration;
        end
    
    end
    
end