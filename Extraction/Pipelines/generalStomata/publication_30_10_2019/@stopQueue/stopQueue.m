classdef stopQueue < eventQueue

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stopQueue constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = stopQueue()
            obj@eventQueue('stopQueue',false);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % do nothing event
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = eventFunction(obj,curNode,eventIdx,dFluid)
            fprintf(['stopQueue event was triggered. Do nothing.\n']);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define trigger - call super class trigger
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = trigger(obj,curNode,dFluid)
            dFluid = trigger@eventQueue(obj,curNode,dFluid);
        end

    end
end