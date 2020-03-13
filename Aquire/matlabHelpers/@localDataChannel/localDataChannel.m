classdef localDataChannel < dataChannel
    properties
        channel;
    end
    
    methods
        
        function [obj] = localDataChannel()
            obj = obj.dataChannel();
            obj.channel = parallel.pool.DataQueue;
        end
        
        function [] = send(obj,data)
            obj.channnel.send(data);
        end
        
    end
    
end