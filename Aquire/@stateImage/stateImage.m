classdef stateImage < handle
    
    properties
       axes;
       rgbStates;
       satLevel;
    end
    
    properties (SetObservable)
        image;
    end
    
    
    methods
        
        function [obj] = stateImage(initImage,rgbStates,satLevel,axes)
            obj.image = initImage;
            obj.axes = axes;
            obj.rgbStates = rgbStates;
            obj.satLevel = satLevel;
        end
        
        function [] = updateTile(obj,r,c,rgb,sat,subImage)
            obj.image = updateStateImage(obj.image,r,c,blockSZ,rgb,sat,subImage,obj.axes);
        end
            
    end
    
end
