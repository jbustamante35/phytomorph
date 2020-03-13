classdef person < doid
    properties
        name;
        email;
    end
    
    methods
        
        function [obj] = person(name,email)
            obj.name = name;
            obj.email = email;
        end
        
    end
  
end