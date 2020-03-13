classdef generalizedFunctionLayer < generalizedDeviceLayer
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor for the generalized function layer
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = generalizedFunctionLayer(name)
            obj = obj@generalizedDeviceLayer(name);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % logger for function latyer
        function [] = logger(obj,dFluid,functionName)
            % super class logger
            logger@generalizedDeviceLayer(obj,functionName);
            % log flow
            dFluid.logger(obj);
        end
        
    end
    
end