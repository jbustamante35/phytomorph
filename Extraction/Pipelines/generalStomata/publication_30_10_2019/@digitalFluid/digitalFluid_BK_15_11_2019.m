classdef digitalFluid < handle & uidTrackable

    properties
        
        % input data moving in the fluid
        data;

        ptr;

        toCompute;
        canParallel;
    


        textHistory;
        logTable;



        flowDirection;

        %%%
        name;
        createDate;
        uid;

        %%%
        persistData;
        stopData;


        jobData;
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % consructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = digitalFluid(name,data)
            

            % human name of fluid
            obj.name = name;
            % fluid create date and time
            obj.createDate = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
            % fluid data
            obj.data = data;



            obj.uid = DataHash([obj.createDate obj.createDate]);
            obj.textHistory = [];
            obj.flowDirection = 'f';

        
            obj.jobData.uid = '1';


            %obj.isTraceRoute = false;
            obj.toCompute = true;

            obj.stopData.uid = '';
            obj.stopData.goFlag = false;

            %{
            if nargin == 2
                obj.stopData.targetQueue = target;
            end
            %}

            obj.persistData.toPersistFunction = @(x)false;
            

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % subref - index the data else normal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'{}')
                varargout = {builtin('subsref',obj.data,S)};
            else
                [varargout{1:nargout}] = builtin('subsref',obj,S);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % subassgn- assign the data else normal
        % if this is not used - should take out
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = subsasgn(obj,S,B)
            if strcmp(S(1).type,'{}')
                obj.data = subsasgn(obj.data,S,B);
            else
                % might have a problem
                obj = builtin('subsasgn',obj,S,B);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the data for the fluid to carry
        % not sure if used - Nov 11 2019
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setData(obj,e,data)
            obj.data{e} = data;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the current
        % not sure if used - Nov 11, 2019
        % looks like the ceLayer
        % may use this method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ret] = getCurrent(obj)
            %ret = obj.data{obj.ptr(1)}{obj.ptr(2)};
            ret = obj.data{obj.ptr(1)};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the current
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ret] = setCurrent(obj)
            %obj.data = subsasgn(obj.data,S,B);
            %ret = obj.data{obj.ptr(1)}{obj.ptr(2)};
            %ret = obj.data{obj.ptr(1)};
        end

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


        function [] = computeToggle(obj,value)
            obj.toCompute = value;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % report the size of the data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [z] = size(obj,varargin)
            z = size(obj.data,varargin{:});
        end

        function [result] = flow(obj,layer)
            % set the target to the layer given to flow
            obj.stopData.targetQueue(1) = layer.bottomCE;
            % call flow
            result = layer.topCE.flow(obj);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % split the stream
        % not sure if used
        % took out Nov 11 2019
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = splitStream(obj,height,width)
            if width > 1
                section = obj.data{height};
                section = repmat({section},[1 width]);
                obj.data{height} = section;
            end
        end


        %{
        function [] = traceRoute(obj,source,target,direction)
            % store state
            
            % suspend computation
            obj.toCompute = false;
            % turn on trace route
            obj.isTraceRoute = true;
            % suspend persistance
            obj.persistData.toPersist = false;

    

            % push new target
            pushTarget(obj,target);
            % flow from source to target
            source.flow(obj);
            % pop target
            popTarget(obj);


            % re-animate computation
            obj.toCompute = true;
            % turn on trace route
            obj.isTraceRoute = false;
            % reanimate the persistance
            obj.persistData.toPersist = false;
        end
        %}


        function [ret] = atCurrentTarget(obj,curLocation)
            ret = obj.stopData.targetQueue(1) == curLocation;
        end

        function [] = pushTarget(obj,newTarget)
            obj.stopData.targetQueue = [newTarget;obj.stopData.targetQueue];
        end

        function [] = popTarget(obj)
            obj.stopData.targetQueue(1) = [];
        end

    end

end