classdef addressBook < oid
    
    properties
        aBook;
    end
    
    methods
        
        function [obj] = addressBook()
   
        end
    
        function [] = add(obj,key,address,channel)
            obj.aBook.(key).address = address;
            obj.aBook.(key).channel = channel;
        end
        
    end
end