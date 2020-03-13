classdef dptr < doid
    properties
        refs;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build out a pointer for the object
        function [obj] = dptr(oin)
            obj = obj@doid();
            if nargin == 1;obj.refs = oid.proj(oin);end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create an oid from the data in refs
        function [refObj] = dereference(obj)
            refObj = oid(obj.refs.type,obj.refs.uuid);
        end
    end
    
    
    methods (Static)
        function [input] = transformInput(input)
            if nargin >= 1
                if ~isempty(input)
                    if ~isa(input,'dptr')
                        input = dptr(input);
                    end
                end
            end
        end
    end
end