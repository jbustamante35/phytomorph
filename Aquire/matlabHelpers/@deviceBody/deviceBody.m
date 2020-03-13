classdef deviceBody < doid
    properties
        
        deviceType;
        resolution;
        timeDelay;
        timeDuration;
        
    end
    
    methods
        function [obj] = deviceBody(deviceType,resolution,timeDelay,timeDuration)
            obj.deviceType = deviceType;
            obj.resolution = resolution;
            obj.timeDelay = timeDelay;
            obj.timeDuration = timeDuration;
        end
    end
end