classdef phytoMessage < doid
    properties
        from;
        to;
        body;
    end
    
    methods
        function [obj] = message(from,to,body)
            if nargin >= 1;obj.from=obj.transformInput(from);end
            if nargin >= 2;obj.to=obj.transformInput(to);end
            if nargin >= 3;obj.body=obj.transformInput(body);end
        end
        
        function [type] = messageType(obj)
            type = obj.body.ref.type;
        end
    end
    
    methods (Access=private)
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