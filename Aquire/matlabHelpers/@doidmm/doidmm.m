classdef doidmm < oidmm
    properties
        genDate;
    end
    
    methods
        function [obj] = doidmm()
            obj = obj@oidmm();
            obj.genDate = datestr(now,'SSMMHHddmmYYYY');
        end
    end
end