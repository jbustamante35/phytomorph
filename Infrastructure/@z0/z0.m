classdef z0 < fa
        
    properties
        
    end
    
    
    methods
    
        function [this] = z0()
            this@fa([]);
        end
        
    end
    
    methods (Static)
        
        function [X] = buildArray(N)
            X = z0.empty(0,N);
            for e = 1:N
                X(e) = z0;
            end
        end
        
    end
end
    