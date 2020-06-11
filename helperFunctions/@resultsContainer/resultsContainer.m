classdef resultsContainer < oid
    % this is the simple object of a set of files.
    % the name generator is special and the index into the fails is special
    % it is a way to hold the output files from many calls toF
    properties (Constant)
        defaultOutputLocation = '/mnt/spaldingdata/nate/octerineDataStore/results/';
    end
    
    properties
        path;
        
        functionName;
        date;
        
        failIndex;
        inputDataCollection;
    end
    
    methods
        function [this] = resultsContainer(projectName,functionName)
            this.date = datestr(now,'yyyy-mm-dd-HH-MM-ss');
            this.functionName = functionName;
            this.path = [resultsContainer.defaultOutputLocation ...
                '{projectName_' projectName '}' filesep ...
                '{functionName_' functionName '}' filesep ...
                '{createDate_' this.date '}' filesep];
            mmkdir(this.path);
        end
        
        function [name] = name(this)
            name = [this.functionName '[' this.inputDataCollection.name ']'];
        end
    end
    
end