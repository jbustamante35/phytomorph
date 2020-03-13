classdef flowable < matlab.mixin.Heterogeneous & handle

    methods (Abstract)
        [result,flowData] = flow(obj,flowData,varargin);
    end

end