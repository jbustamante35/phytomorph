classdef portableQ < doid

    properties
        % bool for open or closed
        isClosed;
        % localQ
        localQ;
        % bool for is connected or free
        isConnected;
        % list of targets the queue feeds into
        targetList;
        % list of sources the queue feeds from
        source;
    end
    
    events
    
    end
    
    methods
    
        function [this] = portableQ(qType)
            
        end
        
        function [] = send(this,data)
            if this.isClosed
                this.localQ{end+1} = data;
            end
        end
        
        function [data,OK] = poll(this)
            if this.isClosed
                if isempty(this.localQ)
                    data = [];
                    OK = false;
                else
                    data = this.localQ{1};
                    data(1) = [];
                    OK = true;
                end
            end
        end
        
        function [] = drain()
            
        end
        
        function [] = QueueLength(this)
            
        end
        
        function [] = addTarget(this,target)
            this.targetList{end+1} = target;
        end
        
        function [] = attachSource(this,source)
            this.source = source;
        end
        
    end

end