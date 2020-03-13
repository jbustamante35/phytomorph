classdef myStateTransitionListener < handle
    
    
    properties
        preTValue;
        postTValue;
        
        preTTarget;
        postTTarget;
        
        passFunc;
    end
    
    properties (SetObservable)
        trigger;
    end
    
  
    methods
        
        function [obj] = myStateTransitionListener(preTTarget,postTTarget,passFunc)
            obj.preTTarget = preTTarget;
            obj.postTTarget = postTTarget;
            obj.passFunc = passFunc;
        end
        
        function [] = addFilterListener(obj,source,property,func)
            % pre capture
            preFunc = @(meta,event)obj.preSet(meta,event);
            % post capture
            postFunc = @(meta,event)obj.postSet(meta,event);
            % attach capture for pre setting
            addlistener(source,property,'PreSet',preFunc);
            % attach capturefor post setting
            addlistener(source,property,'PostSet',postFunc);
        end
        
        % capture the value preset
        function [] = preSet(obj,meta,event)
            obj.preTValue = event.AffectedObject.(meta.Name);
        end
        
        % capture the value post set
        function [] = postSet(obj,meta,event)
            obj.postTValue = event.AffectedObject.(meta.Name);
            preCond = obj.preTValue'*obj.preTTarget;
            postCond = obj.postTValue'*obj.postTTarget;
            if logical(preCond*postCond)
                obj.passFunc(meta,event)
            end
        end
        
        
        
         
        
    end
    
end