classdef functionHandleLayer < generalizedFunctionLayer


    properties
    
        func;

    end


    methods

        %%%%%%%%%%%%%%%%%%%%%%%%
        % constructor for the function handle
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = functionHandleLayer(name,func)
            obj = obj@generalizedFunctionLayer(name);
            obj.func = func;
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%
        % compute the function handle
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [out] = compute(obj,varargin)
            out = obj.func(varargin{:});
            out = {out};
        end

    end

end