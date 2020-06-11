classdef fa < doid
    
    
    properties
        % value of the variable
        value;
        % hash of the data object
        hashValue;
        % file containing value
        dataFile;
    end
    
    methods
    
        
        % constructor
        function [this] = fa(value)
            this.value = value;
            this.hashValue = hash(value);
        end
        
        function [] = save(this)
            global X;
            X = htcArgumentStore();
            this.dataFile = X.add(this);
        end
        
        function [] = clear(this)
            this.value = [];
        end
        
        
        
    end
    
    methods
        
    end
    
end
    