classdef colorConvertLayer < configurableLayer
    
    properties
        type
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%
        % read layer constructor
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = colorConvertLayer(name,type)
            % make layer with name
            obj = obj@configurableLayer(name);
            if nargin == 2;obj.configure(type);end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [] = configure(obj,type)
            obj.type = type;
            %%%%%%%%%%%%%%%%%%%%%%%%
            % set the configuration flag to true
            obj.isConfigured = true;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%
        % compute for the read layer will read the image provided
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [out] = compute(obj,varargin)
            tmp = thawTensor(varargin{3});
            switch obj.type
                case 'rgb>hsv'
                    
                    if size(tmp,3) ~= 3;tmp = cat(3,tmp,tmp,tmp);end

                    tmp = freezeTensor(tmp);
            end
            out = varargin;
            out{4} = tmp;
        end
        
    end
end