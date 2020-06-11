classdef jam < doid
    
    properties (Constant)
        
    end
    
    properties
        
    end
    
    methods (Static)
        
        function [] = spool(o)
            flds = fields(o);
            
            for e = 1:numel(flds)
                flds{e};
            end
        end
    end
    
end