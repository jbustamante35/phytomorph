classdef layerPort < outPort

    
    
    
    methods
        
        function [obj] = layerPort()
            obj@outPort();
            obj.initRemotePath();
        end
        
        function [] = initRemotePath(obj)
            webOptions = weboptions('ContentType','text','Timeout',60);
            %url = 'https://de.cyverse.org/dl/d/E86D57B8-EEF9-4F9F-8CC9-832F04F831DD/networkKey.csv';
            url = 'http://data.cyverse.org/dav-anon/iplant/home/nmiller/publicData/networkKey.csv';
            obj.rPath = webread(url,webOptions);
        end
        
        
    end
    
    methods (Static)
        
        function [] = makeDataPool(poolLocation,poolRoot)
            if nargin == 1;poolRoot = '/iplant/home/nmiller/dataPool/';end
            cmd = ['imkdir -p ' poolRoot poolLocation];
            cmd
            [r,o] = system(cmd);
        end
        
        function [] = removeDataPool(poolLocation,poolRoot)
            if nargin == 1;poolRoot = '/iplant/home/nmiller/dataPool/';end
            cmd = ['irm -r ' poolRoot poolLocation];
            cmd
            [r,o] = system(cmd);
        end
    end
    

end