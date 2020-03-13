classdef sampleMsg < doid
    properties
        file;
        folder;
    end
    
    methods
        
        
        function [obj] = sampleMsg(file,folder)
            obj = obj@doid();
            
            obj.file = file;
            obj.folder = folder;
        end
        
    end
end