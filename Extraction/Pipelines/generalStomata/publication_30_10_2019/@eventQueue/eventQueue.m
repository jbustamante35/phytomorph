classdef eventQueue < matlab.mixin.Heterogeneous  & handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % an eventQueue is a vector of uid(s) 
    % + a vector of event data
    % + an abstract function which is called when one of the uid(s)
    % are equal to the "current" node.uid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties 
        % queue name
        name;
        % node targets
        targetQueue;
        % should clear queue?
        toClear;
        % event data
        eventData;
    end


    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % event queue constructor
        % name := name of event queue
        % toClear := if toClear after trigger
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = eventQueue(name,toClear)
            obj.name = name;
            obj.toClear = toClear;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add node.uid(s) from the queue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = pushEventData(obj,eventData)
            obj.eventData = [eventData obj.eventData];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add node.uid(s) from the queue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = popEventData(obj,ridx)
            if nargin == 1;ridx = 1;end
            obj.eventData(ridx) = [];
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove node.uid(s) from the queue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = pop(obj,ridx)
            if nargin == 1;ridx = numel(obj.targetQueue);end
            obj.targetQueue(ridx) = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add node.uid(s) to the queue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = push(obj,toPush)
            obj.targetQueue{end+1} = toPush;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check the nodes and pop out the true ones
        % if needed to fire
        % note the clearing happens before the firing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = nodeCheck(obj,curNode)
            % for each queue object in the stack
            for e = 1:numel(obj)
                % check the eth node
                b{e} = strcmp(obj(e).targetQueue,curNode.uid);
                % if any nodes should fire and toClear is on for the node
                if (any(b{e}) && any(obj(e).toClear))
                    % find the node to clear
                    if (b{e} && obj.toClear(e))
                        ridx = find(b{e});
                        obj(e).pop(ridx);
                    end
                end
            end
            % if one in stack - unearth the result
            if numel(obj) == 1;b = b{1};end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perform node check and trigger event function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = trigger(obj,curNode,dFluid)
            b = false;
            % for each queue object in the stack
            for e = 1:numel(obj)
                % perform node check
                triggerComplex{e} = obj(e).nodeCheck(curNode);
                % if any are true for "this" node
                b(e) = any(triggerComplex{e});
                if b(e)
                    % find the nodes
                    fidx = find(triggerComplex{e});
                    % run the trigger for each of the nodes which
                    % should fire at this node
                    for n = 1:numel(fidx)
                        dFluid = obj.eventFunction(curNode,fidx(n),dFluid);
                    end
                end
            end
        end


    end

    methods (Abstract)
        
        dFluid = eventFunction(obj,curNode,eventIdx,dFluid)
    
    end

    methods (Sealed)

        function [ret] = find(obj,queueName)
            fidx = strcmp({obj.name},queueName);
            ret = obj(fidx);
        end



    end


end




    

