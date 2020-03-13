classdef initLayer < configurableLayer
    
    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%
        % read layer constructor
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = initLayer(name,configurationData)
            if nargin == 1;configurationData = [];end
            % make layer with name
            obj = obj@configurableLayer(name);
            if ~isempty(configurationData)
                if ischar(configurationData)
                    obj.configure(configurationData);
                else
                    obj.isConfigured = configurationData;
                end
            else
                obj.isConfigured = true;
            end
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%
        % compute for the read layer will read the image provided
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [out] = compute(obj,varargin)
            % userPort = saving the data to the user/condor/test
            % dataPoolPort = data to pour into the community
            % layerPort = layer data for loading and saving
            
            timingBlock('start',['Initializing context@' obj.uid]);
            timingBlock('note','Initializing outPorts {user,data,layer}.');
            global dataPool
            ports = {userPort(varargin{2}),dataPoolPort(),layerPort()};
            portKeys = {'user','data','layer'};
            dataPool = containers.Map(portKeys,ports);
            timingBlock('stop');
            out = varargin;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % configure the read layer with the sample data
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [] = configure(obj,dataSample)
            convertFuntionForm(dataSample);
        end
        
        
        
        
    end
end