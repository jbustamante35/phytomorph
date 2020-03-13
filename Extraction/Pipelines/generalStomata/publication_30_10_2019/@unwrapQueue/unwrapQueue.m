classdef unwrapQueue < eventQueue

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stopQueue constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = unwrapQueue()
            obj@eventQueue('unwrapQueue',false);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % do nothing event
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = eventFunction(obj,curNode,eventIdx,dFluid)
            dFluid = dFluid.data;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define trigger - call super class trigger
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = trigger(obj,curNode,dFluid)
            dFluid = trigger@eventQueue(obj,curNode,dFluid);
        end

    end
end