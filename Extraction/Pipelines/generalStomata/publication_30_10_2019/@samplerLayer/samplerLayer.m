classdef samplerLayer < generalizedFunctionLayer
     
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = samplerLayer(name)
            obj = obj@generalizedFunctionLayer(name);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%
        % sample the image @ the point, over the domain
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [out] = compute(obj,varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%
            toOp = varargin{1};
            pointSet = varargin{2};
            opDomain = varargin{3};
            domainSize = varargin{4};
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % run the sampler
            result = generalizedSampler(toOp,pointSet,opDomain,domainSize);
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % repace the first argument with the sample
            out = varargin;
            out{1} = result;
           
        end
    end
end