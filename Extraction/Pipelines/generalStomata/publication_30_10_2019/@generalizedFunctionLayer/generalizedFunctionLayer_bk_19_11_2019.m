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
        

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pulling the flow call up to the device level - Nov, 11 2019
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow for function layer
        % this should simply call the CE layers above or below depending on
        % the flow direction
        % the compute function will be called if the flow is forward
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = flow(obj,dFluid)
            try
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % flow logic
                %%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(dFluid.flowDirection,'r')


                    % logic of flow through the function layer is passive
                    [dFluid] = obj.topCE.flow(dFluid);

                elseif strcmp(dFluid.flowDirection,'f')


                    % nat can flow through all layers in search of node
                    if dFluid.toCompute
                        % if the forward direction is being called
                        % call compute then flow forwrds via ce
                        %dataIN = dFluid.getCurrent();
                        dataIN = dFluid{1};
                        result = obj.compute(dataIN{:});
                        dFluid{1}{1} = result{1};
                    end


                    % if flow should continue - flow on
                    if ~dFluid.atCurrentTarget(obj)
                        % flow on
                        [dFluid] = obj.bottomCE.flow(dFluid);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%

            catch ME
                report = getReport(ME);
            end
        end
        %}
        
    end
    
end