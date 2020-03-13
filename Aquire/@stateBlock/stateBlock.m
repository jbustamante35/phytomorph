classdef stateBlock < handle &  matlab.mixin.Heterogeneous

    properties
        % determine if/when to move to next block
        transitionFunction;
        % call if transition to next block
        transitionFunctionSuccess;
        % call if fail to transition
        transitionFunctionFail;
        % transformation
        T;
    end
    
    properties (SetObservable)
        % current block state
        currentState;
    end
    
    
    methods
    
        
        function [obj] = stateBlock(T,tF,tfs,tff,initState)
            obj.T = T;
            obj.transitionFunction = tF;
            obj.transitionFunctionSuccess = tfs;
            obj.transitionFunctionFail = tff;
            obj.currentState = initState;
        end
        
        function [] = transition(obj)
            obj.currentState = obj.T*obj.currentState;
        end
        
        function [] = processCommand(obj,meta,event)
            msg = event.AffectedObject.string;
            for e = 1:numel(obj)
                test = obj(e).transitionFunction([],msg);
                vec = obj(e).currentState(:);
                tmp = obj(e).T*kron(test(:),vec);
                tmp = reshape(tmp,[numel(obj(e).currentState) numel(test)]);
                obj(e).currentState = tmp*test';
            end
        end
        
        function [] = setState(obj,state)
            obj.currentState = state;
        end
        
        
    end
    
end