classdef daqMessage < doid & matlab.mixin.SetGet
    properties
        from;
        to;
        body;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%
        % note that the last parameter allows for messages to be used
        % without pointers. This allows for messages to be constructed
        % without the datastore
        %%%%%%%%%%%%%%%%%%
        function [obj] = daqMessage(from,to,body,nonPtr)
            if nargin <=3;nonPtr = false;end
            if ~nonPtr
                if nargin >= 1;obj.from=dptr.transformInput(from);end
                if nargin >= 2;obj.to=dptr.transformInput(to);end
                if nargin >= 3;obj.body=dptr.transformInput(body);end
            else
                if nargin >= 1;obj.from=from;end
                if nargin >= 2;obj.to=to;end
                if nargin >= 3;obj.body=body;end
            end
        end
        
        %%%%%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%%%%
        function [type] = messageType(obj)
            if isa(obj.body,'dptr')
                type = obj.body.refs.type;
            else
                type = obj.body.type;
            end
        end
        
        %%%%%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%%%%
        function [] = set(obj,key,value)
            obj.(key)=dptr.transformInput(value);
        end
    end
    
end