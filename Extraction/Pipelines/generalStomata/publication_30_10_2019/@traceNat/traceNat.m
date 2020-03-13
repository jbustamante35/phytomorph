classdef traceNat < handle

    properties


        isTraceRoute;
        toCompute;

        textHistory;
        flowDirection;

        %%%
        name;
        createDate;
        uid;




        persistData;
        stopData;


        jobData;
    end

    methods


        function [obj] = traceNat(name,target)
            obj.name = name;
            obj.createDate = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
            obj.uid = DataHash([obj.createDate obj.createDate]);
            obj.textHistory = [];
            obj.flowDirection = 'f';

            obj.isTraceRoute = false;
            obj.toCompute = true;

            obj.stopData.uid = '';
            obj.stopData.goFlag = false;
            if nargin == 2
                obj.stopData.targetQueue = target;
            end

            obj.persistData.toPersist = false;

            obj.jobData.uid = '1';

        end

        % log route
        function [] = logLocation(obj,curLocation)
            logString = ...
            string(['{class_' class(curLocation) ...
                '}{name_' curLocation.name ...
                '}{uid_' curLocation.uid ...
                '}{date_' datestr(datetime,'ss_hh_dd_mm_yyyy') ...
                '}{direction_' obj.flowDirection '}']);
            obj.textHistory = [obj.textHistory;logString];
        end




        function [] = pushState(obj,stateVector)

        end
        

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

        function [ret] = atCurrentTarget(obj,curLocation)
            ret = obj.stopData.targetQueue(1) == curLocation;
        end

        function [] = pushTarget(obj,newTarget)
            obj.stopData.targetQueue = [newTarget;obj.stopData.targetQueue];
        end

        function [] = popTarget(obj)
            obj.stopData.targetQueue(1) = [];;
        end

    end

end