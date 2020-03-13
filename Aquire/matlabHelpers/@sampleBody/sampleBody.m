classdef sampleBody < doid
    properties
        file;
        folder;
    end
    
    methods
        
        
        function [obj] = sampleBody(file,folder)
            obj = obj@doid();
            
            obj.file = file;
            obj.folder = folder;
        end
        
    end
end