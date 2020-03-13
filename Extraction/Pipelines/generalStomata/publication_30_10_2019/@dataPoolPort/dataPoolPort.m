classdef dataPoolPort < outPort

    
    
    
    methods
        
        function [obj] = dataPoolPort()
            obj@outPort();
            obj.initRemotePath();
        end
        
        
        %
        function [] = initRemotePath(obj)
            webOptions = weboptions('ContentType','text','Timeout',60);
            %url = 'https://de.cyverse.org/dl/d/88EE12DA-A675-4234-93F4-48B6E3E147EA/poolKey.csv';
            url = 'http://data.cyverse.org/dav-anon/iplant/home/nmiller/publicData/poolKey.csv';
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