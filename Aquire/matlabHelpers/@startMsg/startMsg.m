classdef startMsg < doid
    properties
        deviceType;
        resolution;
        timeDelay;
        timeDuration;
        
        
        file;
        folder;
    end
    
    methods
        
        
        function [obj] = startMsg(file,folder,deviceType,resolution,timeDelay,timeDuration)
            obj = obj@doid();
            
            
            obj.file = file;
            obj.folder = folder;
            
            
            obj.deviceType = deviceType;
            obj.resolution = resolution;
            obj.timeDelay = timeDelay;
            obj.timeDuration = timeDuration;
        end
        
    end
end