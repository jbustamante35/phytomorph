classdef addressBook < doid
    
    properties
        entries;
    end
    
    
    methods
        % constructor for address book
        % init with self
        function [obj] = addressBook(self_uuid,self_channel)
            obj = obj@doid();
            selfEntry = addressBookEntry(self_uuid,self_channel,{});
            obj.addEntry(selfEntry);
        end
        
        % add entry
        function [] = addEntry(obj,entries)
            obj.entries = [obj.entries;entries];
        end
        
        function [entry] = getSelfEntry(obj)
            entry = obj.entries(1);
        end
        
        
        function [] = isCollision(obj,entry)
            here = 1;
        end
        
        
    end
    
end