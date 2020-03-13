classdef digitalFluid < handle & uidTrackable & inputBlock

    properties
    

        %%%%%%%%%%%%%%%%%%%%%%%%
        % compute
        toCompute;
        canParallel;


        %%%%%%%%%%%%%%%%%%%%%%%%
        % text log of flow through network
        textHistory;

        %%%%%%%%%%%%%%%%%%%%%%%%
        % direction and color of flow for the fluid
        flowDirection;
        color;


        %%%%%%%%%%%%%%%
        % stackable queueManifold for event triggers
        queueManifold;
    

        jobData;

        errorLog;

    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % consructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = digitalFluid(name,data,color)
            
            obj = obj@uidTrackable(name);
            obj = obj@inputBlock(data);

            if nargin ~=3;color = "green";end


            obj.textHistory = [];
            obj.flowDirection = 'f';
            obj.color = color;

        
            obj.jobData.uid = '1';


            obj.toCompute = true;


            obj.queueManifold = [persistQueue();stopQueue();unwrapQueue()];

        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the current
        % not sure if used - Nov 11, 2019
        % looks like the ceLayer
        % may use this method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [data] = getCurrent(obj)
            data = obj.data;
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the current
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setCurrent(obj,data)
            obj.data = data;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trigger all/one/set of event queue(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = trigger(obj,curNode,eventIndex)
            if nargin == 2;eventIndex = 1:numel(obj.queueManifold);end
            for e = 1:numel(eventIndex)
                obj = obj.queueManifold(eventIndex(e)).trigger(curNode,obj);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add persist location
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = addPersistNode(obj,pNode,eventData)
            persistNode = obj.queueManifold.find('persistQueue');
            persistNode.pushEventData(eventData);
            persistNode.push(pNode);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % is all fluid(s) are green - true inputs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isGreen(obj)
             b = obj.isColor("green");
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % is all fluid(s) are blue - complex of (r,g,b)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isBlue(obj)
             b = obj.isColor("blue");
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % is all fluid(s) are red - persisted data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isRed(obj)
            b = obj.isColor("red");
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % is all fluid(s) red wrapped in blue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isPurple(obj)
            b = obj.isBlue();
            if isa(obj.data,'digitalFluid')
                b = b & obj.data.isRed();
            else
                b = false;
            end
        end


        function [] = red2green(obj)
            tmpFluid = generalizeLoader(obj.data);
            obj = tmpFluid.dFluid;
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % general color match
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isColor(obj,colorMatch)
            b = all(strcmp([obj.color],colorMatch));
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % log fluidflow through pipeline
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = logger(obj,curLocation)
            logString = ...
            string(['{class_' class(curLocation) ...
                '}{name_' curLocation.name ...
                '}{uid_' curLocation.uid ...
                '}{date_' datestr(datetime,'ss_hh_dd_mm_yyyy') ...
                '}{direction_' obj.flowDirection '}']);
            obj.textHistory = [obj.textHistory;logString];
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % report the size of the data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = computeToggle(obj,value)
            obj.toCompute = value;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if fluid is wide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = isWide(obj)
            b = size(obj,2) > 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow short-cut
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [result] = flow(obj,layer,direction)
            % set default values for direction as forward 'f'
            if nargin == 2;direction = 'f';end
            if strcmp(direction,'r');obj.flowDirection = direction;end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if calling on a 'true' layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isa(layer,'ceLayer')
                % push the stop location to the manifold of queues
                obj.queueManifold(2).push(layer.bottomCE.uid);
                % call flow
                result = layer.topCE.flow(obj);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calling with direct instructions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                % queue the given id to the stop(2) queue
                obj.queueManifold(2).push(layer.uid);
                % call flow
                result = layer.flow(obj);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syntax sugar for stopping at a target node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [b] = toStop(obj,curNode)
             b = obj.queueManifold(2).nodeCheck(curNode);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syntax sugar for adding a stop target
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = pushStopTarget(obj,newTarget)
            obj.queueManifold(2).push(newTarget.uid);
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syntax sugar for removing a stop target
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = popStopTarget(obj)
            obj.queueManifold(2).pop();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syntax sugar for adding an unwrap point for blue fluid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = assignUnwrapNode(obj,unwrapTarget)
            obj.queueManifold(3).push(unwrapTarget.uid);
        end
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % numel of targets in queue size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [s] = queueSize(obj)
            for e = 1:numel(obj.queueManifold)
                s(e) = numel(obj.queueManifold(e).targetQueue);
            end
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

end