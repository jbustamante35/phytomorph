classdef generalizedDeviceLayer < generalizedComputeLayer & linkable

    properties
        topCE;
        bottomCE;



        localGraph;
        nodeList;
        connectionTensor;
    end
    
    methods
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor the the device layer
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = generalizedDeviceLayer(name)
            % call super constructur
            obj = obj@generalizedComputeLayer(name);
            % call the linkable super constructor
            obj = obj@linkable();
            % construct the top collector/emitter layer
            obj.topCE = ceLayer([name '-topCE'],'top');
            % construct the bottom collector/emitter layer
            obj.bottomCE = ceLayer([name '-bottomCE'],'bottom');
            % link-aware the top to the bottom
            obj.topCE.setPairLink(obj.bottomCE);
            % link-aware the bottom to the top
            obj.bottomCE.setPairLink(obj.topCE);


            % link the device layer upwards to the topCE
            obj.linkUp(obj.topCE);
            % link the topCE to the device layer
            obj.topCE.linkDown(obj);
            % link the device layer downwards to the bottomCE
            obj.linkDown(obj.bottomCE);
            % link the bottomCE upward to the device layer
            obj.bottomCE.linkUp(obj);




            % update the graph for this object
            obj.updateGraph();
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % not needed here?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function result = compute(obj,varargin)
            result = {};
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call view of the graph-map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = view(obj)
            view(obj.localGraph);
        end

        function [returnSize] = graphSize(obj)
            sz = [3 1 1 1];
            returnSize.sz = sz;
            returnSize.partition = [1];
        end

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
                    if (dFluid.toCompute && dFluid.isGreen())
                        % if the forward direction is being called
                        % call compute then flow forwards via ce
                        
        
                        dataIN = dFluid.getCurrent();
                        %dataIN = dFluid{1};
                        

                        dataOUT = obj.compute(dataIN{:});
                        
                        dFluid.setCurrent(dataOUT);
                        %dFluid{dFluid.ptr} = result;
                        %dFluid{1}{1} = result{1};
                    end


                    % if flow should continue - flow on
                    %if ~dFluid.atCurrentStopTarget(obj)
                    if ~dFluid.toStop(obj)
                        % flow on
                        [dFluid] = obj.bottomCE.flow(dFluid);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%

            catch ME
                dFluid.errorLog = getReport(ME);
            end
        end
        
        %{
        function [] = viewW(obj,cons)
            
            

            % n1 = obj.NlocalGraph.getnodesbyid(obj.NnodeList{1});

            % bg.nodes(1).Positio

            %sz = graphSize(obj);
    
        end
        %}

        function [] = updateGraph(obj)


            [nodeList,upid,M] = obj.topCE.getNodeList();
            upid(1,:) = [];
            connectionMatrix = zeros(numel(nodeList));

            for e = 1:size(upid,1)
                childIDX = strcmp(nodeList,upid{e,1});
                parentIDX = strcmp(nodeList,upid{e,2});
                connectionMatrix(parentIDX,childIDX) = 1;
            end


            obj.localGraph = biograph(connectionMatrix,nodeList);
            %obj.nodeManifest = M;
            obj.nodeList = nodeList;


            for e = 1:numel(nodeList)
                if contains(nodeList{e},'ceLayer')
                    obj.localGraph.nodes(e).Shape = 'circle';
                    obj.localGraph.nodes(e).Size = [2 2];
                    obj.localGraph.nodes(e).Label = 'ceLayer';
                else
                    obj.localGraph.nodes(e).Label = nodeList{e};
                end
            end


            %{
            fluxMatrix = [[0 1 0];[0 0 1];[0 0 0]];
    
            nodeList{1} = obj.topCE.getUidType();
            nodeList{2} = obj.getUidType();
            nodeList{3} = obj.bottomCE.getUidType();



            %obj.NlocalGraph = biograph(fluxMatrix,nodeList,'CustomNodeDrawFcn',@(X)generalizedDeviceLayer.view(X));
            obj.localGraph = biograph(fluxMatrix,nodeList);
         

            n1 = obj.localGraph.getnodesbyid(nodeList{1});
            n2 = obj.localGraph.getnodesbyid(nodeList{2});
            n3 = obj.localGraph.getnodesbyid(nodeList{3});
            
            n1.Label = 'ceLayer';
            n2.Label = obj.getUidType();
            n3.Label = 'ceLayer';

            n1.Shape = 'circle';
            n3.Shape = 'circle';
          
            
            

            obj.nodeList = nodeList;
            obj.connectionTensor = fluxMatrix;
            %}
        end

        
    end

    %{
    methods (Static)

        function [Y] = view(X)
        
         here = 1;

        end
    end
    %}
end