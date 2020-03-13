classdef msgShredder < handle

    
    properties
    
        toFunc;
        fromFunc;
        msgFunc;
        
    end
    
    
    methods
        function [obj] = msgShredder(toFunc,fromFunc,msgFunc)
            obj.toFunc = toFunc;
            obj.fromFunc = fromFunc;
            obj.msgFunc = msgFunc;
        end
        
        function [] = shred(obj,msg)
            
        end
        
    end
    
    
    
    
end