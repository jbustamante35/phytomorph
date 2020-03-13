classdef loopStack < generalizedDeviceLayer

    properties
        layerStack;
        totalLoops = N;
    end

    methods
        function [obj] = loopStack(name,layerStack)
            % call super constructur
            obj = obj@generalizedDeviceLayer(name);


            obj.layerStack = layerStack;

            % connect the CEtop to this - compute/flow
            obj.topCE.downList{1} = obj.layerStack.topCE;


            % connect the CEbottom to this - compute/flow
            obj.bottomCE.upList{1} = obj.layerStack.topCE;


        end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow for function layer
        % this should simply call the CE layers above or below depending on
        % the flow direction
        % the compute function will be called if the flow is forward
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [result,flowData] = flow(obj,flowData,varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % log flow
            flowData.logLocation(obj);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % flow logic
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(flowData.flowDirection,'r')

        
                % logic of flow through the function layer is passive
                [result,flowData] = obj.topCE.flow(flowData,varargin{:});


            elseif strcmp(flowData.flowDirection,'f')


                % nat can flow through all layers in search of node
                if flowData.toCompute
                    % if the forward direction is being called
                    % call compute then flow forwrds via ce
                    result = obj.compute(varargin{:});
                else
                    % do not compute
                    result = {};
                end


                % if flow should continue - flow on
                if ~flowData.atCurrentTarget(obj)
                    % flow on
                    [result,flowData] = obj.bottomCE.flow(flowData,result{:});
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%


        end
        

        

    end

end