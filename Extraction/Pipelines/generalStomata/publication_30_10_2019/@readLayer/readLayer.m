classdef readLayer < configurableLayer
    
    properties
        readParameters;
        type;
        cachedValue;
        cachedKey;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%
        % read layer constructor
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = readLayer(name)
            % make layer with name
            obj = obj@configurableLayer(name);
            % set configure to false
            obj.isConfigured = false;
            
            % init cached data;
            obj.cachedValue = [];
            % init cached name
            obj.cachedKey = '';
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%
        % compute for the read layer will read the image provided
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [out] = compute(obj,varargin)
            timingBlock('start',['Image read layer@:' obj.uid]);
            % first input must be the file name
            fileName = varargin{1};
            % if the object configured then continure
            if obj.isConfigured
                % if the requested image is cached
                if strcmp(fileName,obj.cachedKey)
                    timingBlock('note','Loading cached copy.');
                    % return the cached copy
                    result = obj.cachedValue;
                else
                    timingBlock('note','Loading disk copy.');
                    if ~isempty(obj.readParameters)
                        % if not cached then read and cache the image
                        [result] = generalizeLoader(fileName,obj.type,obj.readParameters(1));
                    else
                        [result] = generalizeLoader(fileName,obj.type);
                    end
                    % freeze back to format for pasing to next function
                    result = freezeTensor(result);
                    % set the cached copy
                    obj.cachedKey = fileName;
                    obj.cachedValue = result;
                end
                % set the output and return
                out = varargin;
                out{3} = result;
            else
                % attempt to auto configure with the dataset
                fprintf(['Please provide data to configure.\n']);
            end
            timingBlock('stop');
        end
        
        function [] = setType(obj,type)
            obj.type = type;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % configure the read layer with the sample data
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [] = configure(obj,dataSample)
            if ~obj.isConfigured
                %%%%%%%%%%%%%%%%%%%%%%%%
                % extract only the image names from the datasample
                %%%%%%%%%%%%%%%%%%%%%%%%
                for s = 1:numel(dataSample)
                    imageNames{s} = dataSample{s}{1};
                end
                %%%%%%%%%%%%%%%%%%%%%%%%
                % configure the read layer
                obj.readParameters{1} = gatherParametersForReader(imageNames,obj.type,false);
                %%%%%%%%%%%%%%%%%%%%%%%%
                % set the configuration flag to true
                obj.isConfigured = true;
            end
        end
        
    end
end