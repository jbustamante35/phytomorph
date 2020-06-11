classdef tracer < doidmm
    properties
        data;
    end
    
    methods
        
        function [this] = tracer(data)
            this.data = data;
        end
        
        function [c] = class(this)
            c = class(this.data);
        end
        
        function [varargout] = subsref(this,subs)
            if nargout == 0
                builtin('subsref',this.data,subs)
            else
                varargout{1:nargout} = builtin('subsref',this.data,subs);
            end
        end
        
        
        
        
    end
    
end