classdef resultsCollection < dataCollection
    
    properties (Constant)
        outportLocation = '/mnt/spaldingdata/nate/octerineDataStore/results/outport/';
    end

    properties
        fsStatus;
    end
    
    methods
        
        function [obj] = resultsCollection(fileList)
            obj@dataCollection(fileList);
        end
        
        function [] = setStatus(this,status)
            this.fsStatus = status;
        end
        
        function [] = copyResults(this,target)
            for e = 1:numel(this)
                for f = 1:numel(this(e).FileList)
                    this(e).FileList{f}.copyFile(target);
                end
            end
        end
        
        function [fileList] = findImageFiles(this)
            fileList = {};
            for e = 1:numel(this)
                fileList{e} = {};
                for i = 1:numel(this(e).FileList)
                    if this(e).FileList{i}.isImage
                        fileList{e}{end+1} = this(e).FileList{i};
                    end
                end
            end
            if numel(this) == 1
                fileList = fileList{1};
            end
        end
        
    end
    
    methods (Static)
        
        function [splitLocation] = splitResults(baseLocation,zipName)
            % init the data collection
            resultCollection = dataCollection.buildFromDirectory(baseLocation);
            % make a temp location
            tLoc = tempname;
            tLoc = tLoc(6:end);
            tLoc = [resultsCollection.outportLocation tLoc];
            zipFile = [tLoc filesep zipName '.zip'];
            % make a success location under temp
            succussLoc = [tLoc filesep 'success' filesep];
            % make a fail location under temp
            failLoc = [tLoc filesep 'fail' filesep];
            mmkdir(failLoc);
            mmkdir(succussLoc);
            % find the fail
            fidx0 = strcmp({resultCollection.fsStatus},'fail');
            % find the success
            fidx1 = strcmp({resultCollection.fsStatus},'success');
            % make a results collection of fails
            resultCollection0 = resultCollection(fidx0);
            % make a results collection of successes
            resultCollection1 = resultCollection(fidx1);
            % copy the fails
            resultCollection0.copyResults(failLoc);
            % copy the successes
            resultCollection1.copyResults(succussLoc);
            % find zip files
            zipList = dig(tLoc);
            % zip up data
            zip(zipFile,zipList);
        end
        
        function [] = view(baseLocation,type)
            if nargin == 0;type = 'images';end
            % init the data collection
            resultCollection = dataCollection.buildFromDirectory(baseLocation);
            % loop over each collection
            for e = 1:numel(resultCollection)
                imageFiles = resultCollection(e).findImageFiles();
                for i = 1:numel(imageFiles)
                    I = imread(imageFiles{i});
                    imshow(I,[]);
                    drawnow
                    waitforbuttonpress
                end
                %[d,A]=strdist(r,b,krk,cas)
            end
        end
        
    end
    
end