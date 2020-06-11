classdef htcDag < doid

    
    properties
        nodeList;
        edgeList;
        
    end
    
    
    methods
        function [this] = htcDag()
            
        end
        
        function [] = addNode(this,node)
            this.nodeList{end+1} = node.uuid;
        end
        
        function [] = addEdge(this,source,target)
            sourceIdx = this.search(source);
            targetIdx = this.search(target);
            this.edgeList{end+1} = [sourceIdx targetIdx];
        end
        
        function [idx] = search(this,node,toInsert)
            if nargin == 2;toInsert = true;end
            idx = find(strcmp(this.nodeList,node.uuid));
            if isempty(idx) && toInsert
                this.addNode(node);
                idx = this.numelNodes();
            end
        end
        
        function [n] = numelNodes(this)
            n = numel(this.nodeList);
        end
        
        function [] = generateDAGfile(this)
            here = 1;
        end
        
    end
    
end