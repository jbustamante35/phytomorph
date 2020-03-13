classdef layerStack < configurableLayer
    
    
    properties
        layers;
        layersMap;
        nodeManifest;
        nodeList;
    end
    
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = layerStack(name,layers,glueCommands)
            
            obj = obj@configurableLayer(name);
            obj.layers = layers;


            obj.topCE.downList = {};
            obj.topCE.upList = {};
           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if needed - hook up funnel layer to this.bottomCE
            if (size(obj.layers,2) > 1)

                %if nargin == 2
                %    glueCommands = ceLayer.defaultCommands(size(obj.layers,2));
                %end

                %obj.bottomCE.setFunnelCommands(glueCommands);

                % create funnel layer
                %dataFunnel = funnelLayer([name '--dataFunnel'],glueCommands);



                 %obj.bottomCE.upList{1} = dataFunnel.bottomCE;
                 %dataFunnel.bottomCE.downList{1} = obj.bottomCE;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



            for w = 1:size(obj.layers,2)

                size(obj.layers,2) ;
                if size(obj.layers,2) > 1
                    here = 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the list of CE nodes to hook up
                bottomTMPlist = [obj.topCE,obj.layers(:,w).bottomCE];
                topTMPlist = [obj.layers(:,w).topCE,obj.bottomCE];
                % note that this will wire different than the above two
                %bottomTMPlist = [obj.layers(1:end-1,w).bottomCE];
                %topTMPlist = [obj.layers(2:end,w).topCE];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % hook up the required nodes
                for e = 1:numel(topTMPlist)
                    % hook the bottom (starting with this.topCE) to top
                    %bottomTMPlist(e).downList{w} = topTMPlist(e);
                    bottomTMPlist(e).downList{end+1} = topTMPlist(e);

                    %{
                    if isempty(bottomTMPlist(e).downList{w})
                        here = 1;
                    end
                    %}

                    topTMPlist(e).upList{end+1} = bottomTMPlist(e);

                    %{
                    if e == numel(topTMPlist)
                        % reverse
                        topTMPlist(e).upList{w} = bottomTMPlist(e);
                    else
                        % reverse
                        topTMPlist(e).upList{1} = bottomTMPlist(e);
                    end
                    %}
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

                if (size(obj.layers,2) > 1)


                   
                    %{
                    % removed Nov 11, 2019 - put back Nov, 11 2019
                    % removed again Nov 11, 2019 for different reason
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % upper CE to the top CE of the eth layer in the stack
                    obj.topCE.downList{w} = obj.layers(1,w).topCE;
                    % reverse the connection
                    obj.layers(1,w).topCE.upList{1} = obj.topCE;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %}


                    %{
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % bottom CE to the bottom of each in the stack
                    obj.bottomCE.upList{w} = obj.layers(end,w).bottomCE;
                    % reverse the connection
                    obj.layers(end,w).bottomCE.downList{1} = obj.bottomCE;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %}

                    
                    %{
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % bottom CE to the bottom of each in the stack
                    dataFunnel.topCE.upList{w} = obj.layers(end,w).bottomCE;
                    % reverse the connection
                    obj.layers(end,w).bottomCE.downList{1} = dataFunnel.topCE;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %}
                else
                    %{
                    % removed Nov 11, 2019 - put back Nov, 11 2019
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % upper CE to the top CE of the eth layer in the stack
                    obj.topCE.downList{w} = obj.layers(1,w).topCE;
                    % reverse the connection
                    obj.layers(1,w).topCE.upList{1} = obj.topCE;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % bottom CE to the bottom of each in the stack
                    obj.bottomCE.upList{w} = obj.layers(end,w).bottomCE;
                    % reverse the connection
                    obj.layers(end,w).bottomCE.downList{1} = obj.bottomCE;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %}
                end


            end





            %{
            if size(obj.layers,1) > 1

                %{
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the list of CE nodes to hook up
                bottomTMPlist = [obj.topCE,obj.layers.bottomCE];
                topTMPlist = [obj.layers.topCE,obj.bottomCE];


                 % hook up the required nodes
                for e = 1:numel(topTMPlist)
                    bottomTMPlist(e).downList{1} = topTMPlist(e);
                    topTMPlist(e).upList{1} = bottomTMPlist(e);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % wide stack
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





            elseif size(obj.layers,2) > 1
               

                %{
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % attach each wide layer ce to this ce top
                for e = 1:size(obj.layers,2)
                    % upper CE to the top CE of the eth layer in the stack
                    obj.topCE.downList{e} = obj.layers(e).topCE;
                    % reverse the connection
                    obj.layers(e).topCE.upList{1} = obj.topCE;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % attach each wide layer ce to this ce bottom
                for e = 1:size(obj.layers,2)
                    % bottom CE to the bottom of each in the stack
                    obj.bottomCE.upList{e} = obj.layers(e).bottomCE;
                    % reverse the connection
                    obj.layers(e).bottomCE.downList{1} = obj.bottomCE;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}
            end


            %}

            obj.updateLayersMap();
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % call configure with the sample data to each of the layers
        %%%%%%%%%%%%%%%%%%%%
        function [] = configure(obj,dataSample)
            for e = 1:numel(obj.layers)
                if isa(obj.layers(e),'configurableLayer')
                    obj.layers(e).configure(dataSample);
                end
            end
        end
        

        function [] = getDependancyNodeList(obj,sourceNode,targetNode)
            if nargin == 2
                targetNode = obj.topCE.getUidType();
                targetNode = obj.NodeToNum(targetNode);
            end
            sourceNode = obj.NodeToNum(sourceNode);
            [pDist,nodePath] = obj.layersMap.shortestpath(targetNode,sourceNode);
        end



        function [] = updateLayersMap(obj)

            fluid = digitalFluid('mapFluid',{});
            fluid.computeToggle(false);
            res1 = fluid.flow(obj);




            [nodeList,upid,M] = obj.topCE.getNodeList();
            upid(1,:) = [];
            connectionMatrix = zeros(numel(nodeList));

            for e = 1:size(upid,1)
                childIDX = strcmp(nodeList,upid{e,1});
                parentIDX = strcmp(nodeList,upid{e,2});
                connectionMatrix(parentIDX,childIDX) = 1;
            end


            obj.layersMap = biograph(connectionMatrix,nodeList);
            obj.nodeManifest = M;
            obj.nodeList = nodeList;

            for e = 1:numel(nodeList)
                if contains(nodeList{e},'ceLayer')
                    obj.layersMap.nodes(e).Shape = 'circle';
                    obj.layersMap.nodes(e).Size = [2 2];
                    obj.layersMap.nodes(e).Label = 'ceLayer';
                else
                    obj.layersMap.nodes(e).Label = nodeList{e};
                end
            end

        end

        function [] = view(obj)

            view(obj.layersMap);

        end


        
        function [idx] = NodeToNum(obj,layerUidType)
            idx = find(strcmp(obj.nodeList,layerUidType));
        end
    

        
        function [result] = compute(obj,varargin)
            % As of Nov 12, 2019
            % this layer should not need to compute
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if this layer is asked to compute then
            % as a stack layer it will call flow over the layer
            % stack
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %{
            % get the height of the data
            inputStackHeight = size(dFluid,1);
            % get the width of the data
            Width = numel(obj.downList);

            wFluid = digitalFluid('workingFluid',{varargin});
            wFluid.stopData.targetQueue = dFluid.stopData.targetQueue;
            wFluid.toCompute = dFluid.toCompute;

            %}


        end
    end
    

    
end