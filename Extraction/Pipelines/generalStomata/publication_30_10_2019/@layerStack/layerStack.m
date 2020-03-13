classdef layerStack < configurableLayer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % layerStack object := defaults to a configurable layer
    % as a connfigurable-layer, these are stackable devices.
    % These are the first complex building block above a simple device
    % calling the constructor with an array of layers will start the
    % gluing process. Therefore anytime layers need to be stacked and
    % glued, this is the constructor to call.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % constructor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % the input matrix of layers
        layers;
        % key-value map for the layers in the matrix
        layersMap;
        
        nodeManifest;


    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        % name := human name of layer
        % layers := matrix of layers to glue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = layerStack(name,layers,glueCommands)
            % call the super constructor
            obj = obj@configurableLayer(name);
            % set the layers matrix
            obj.layers = layers;
            % clear links to the layer
            obj.clearAllLinks();
            obj.topCE.clearAllLinks();
            obj.bottomCE.clearAllLinks();
            % loop over the width of the input stack
            for w = 1:size(obj.layers,2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the list of CE nodes to hook up
                % array of bottom nodes - first bottom node is this.topCE
                % this is bottom from layer aboves perspective
                bottomTMPlist = [obj.topCE,obj.layers(:,w).bottomCE];
                % array of top nodes to hook to the bottom
                topTMPlist = [obj.layers(:,w).topCE,obj.bottomCE];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % hook up the required nodes
                for e = 1:numel(topTMPlist)
                    % add the top element from the node below to the stack
                    bottomTMPlist(e).linkDown(topTMPlist(e));
                    % link the bottom element to the uplist
                    topTMPlist(e).linkUp(bottomTMPlist(e));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            % update the layers map

            % update the graph for this object
            obj.updateGraph();
            %obj.updateLayersMap();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


        
        function [returnSize] = graphSize(obj)
            % for each column in the layer matrix
            for c = 1:size(obj.layers,2)
                % init the width for each column to zero
                colW(c) = 0;
                colH(c) = 0;
                % loop down the row of the c-th column
                for r = 1:size(obj.layers,1)
                    % column is the max width of the r-th row
                    tmpSZ{r,c} = obj.layers(r,c).graphSize;
                    colH(c) = colH(c) + tmpSZ{r,c}.sz(1);
                    colW(c) = max(colW(c),tmpSZ{r,c}.sz(2));
                end
                % add extra for ceBlocks at the top/bottom of layerStack
                colH(c) = colH(c) + 2;
            end
            % size and blocking of layers
            sz = [max(colH) sum(colW) size(obj.layers)];
            returnSize.sz = sz;
            returnSize.partition = tmpSZ;
        end

        function [] = viewW(obj,extra)
            perUnit = 100;
            padUnit = .05*perUnit;
            sz = obj.graphSize();



            if nargin == 1;extra.fig = figure;extra.const = [0 0 sz.sz(2)*perUnit sz.sz(1)*perUnit];end
            





            rectangle('Position',extra.const,'FaceColor','g');
            hold all
            axis equal
            
        
            extra.const(1) = extra.const(1) + padUnit;
            

            for r = 1:size(obj.layers,1)
        
                %extra.const(2) = extra.const(2) + perUnit*1;
                %extra.const(4) = extra.const(4) - perUnit*2;

                for c = 1:size(obj.layers,2)
                    %{
                    extra.const(1) = extra.const(1) + perUnit*(c-1);

                    extra.const(3) = extra.const(3) - padUnit*2;

                    obj.layers(r,c).viewW(extra);
                    %}
                    tmpSZ = obj.layers(r,c).graphSize();


                    
                    
                end
            end


            %{
            for e = 1:numel(obj.layers)
            end
            %}
            here =1;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call configure with the sample data to each of the layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = configure(obj,dataSample)
            for e = 1:numel(obj.layers)
                if isa(obj.layers(e),'configurableLayer')
                    obj.layers(e).configure(dataSample);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        function [] = getDependancyNodeList(obj,sourceNode,targetNode)
            if nargin == 2
                targetNode = obj.topCE.getUidType();
                targetNode = obj.NodeToNum(targetNode);
            end
            sourceNode = obj.NodeToNum(sourceNode);
            [pDist,nodePath] = obj.layersMap.shortestpath(targetNode,sourceNode);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % walk the layers stack and generate the
        function [] = updateLayersMap(obj)


            % return the node list, uniqueID->parentID cell array and
            % 
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


        

      
    

        
        function [result] = compute(obj,varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % As of documentation date Nov 15, 2019
            % this layer should not need to compute
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        end
    end
    

    methods (Access = private)
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert a layer name to its number name in the list of nodes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [idx] = NodeToNum(obj,layerUidType)
            idx = find(strcmp(obj.nodeList,layerUidType));
        end

    end

    
end