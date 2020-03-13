classdef fftLayer < generalizedFunctionLayer
    
    properties
        frequencyToExtract;
    end
    
    methods
        function [obj] = fftLayer(name)
            obj = obj@generalizedFunctionLayer(name);
        end
        
        function [] = setN(obj,N)
            obj.frequencyToExtract = N;
        end
        
        
        function [out] = compute(obj,varargin)
            patchData = thawTensor(varargin{1},2);
            M = mfftm(size(patchData,2),obj.frequencyToExtract);
            patchData = permute(patchData,[2 1]);
            fF = mtimesx(M,patchData);
            fF = permute(fF,[2 1]);
            
            T{1} = abs(fF);
            T{2} = angle(fF);
            

            out = varargin;
            subOut{1} = patchData;
            subOut{2} = fF;
            subOut{3} = abs(fF);
            subOut{4} = angle(fF);
            subOut = freezeTensor(subOut);

            out{1} = subOut;

        end
    end
end