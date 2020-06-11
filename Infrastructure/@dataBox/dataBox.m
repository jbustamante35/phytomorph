classdef dataBox
    
    properties
        dataLocations;
    end
    
    methods
        
        function [obj] = dataBox(dataLocations)
            if iscell(dataLocations)
                obj.dataLocations = dataLocations;
            else
                obj.dataLocations = {dataLocations};
            end
        end
        
    end
    
end