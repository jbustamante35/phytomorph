classdef doid < oid
    properties
        genDate;
    end
    
    
    methods
        function [obj] = doid()
            obj = obj@oid();
            obj.genDate = datestr(now,'SSMMHHddmmYYYY');
        end
    end
end