classdef testable < handle

    properties
        testName;
    end
    
    methods 
        function [obj] = testable(name)
            obj.testName = name;
        end
    end

end