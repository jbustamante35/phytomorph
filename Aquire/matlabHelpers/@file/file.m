classdef file < doid

    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % path to file
        filePath;
        % name of file
        fileName;
        % file size
        size;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % hash of data
        dataHash;
        % hash of name
        nameHash;
        % hash of data+name
        fullHash;
        % might now be needed
        stage = 1;
    end

    methods
        function [obj] = file(fName,metaData)
            % super constructor
            obj = obj@doid();
            if nargin ~= 0
                
                % get the file parts
                [obj.filePath,nm,ext] = fileparts(fName);
                obj.filePath = [obj.filePath filesep];
                obj.fileName = [nm ext];



                if nargin == 1
                    [obj.size,obj.dataHash] = obj.hashData();
                else
                    obj.size = metaData.size;
                    obj.dataHash = metaData.dataHash;
                end
                obj.nameHash = obj.hashName();
                obj.fullHash = obj.hashAll();


                % overwrite the uuid as the full-hash=hash(hash(data)+hash(name))
                obj.uuid = obj.fullHash;
            end
        end
        
        function [sz,hash] = hashData(obj)
            tmpName = [obj.filePath obj.fileName];
            if isLOCAL(tmpName)
                [sz,hash] = file.getLocalproperties(tmpName);
            elseif isCHTC(tmpName)
                [sz,hash] = file.getCHTCproperties(tmpName);
            elseif isIRODS(tmpName)
                [sz,hash] = file.getIRODSproperties(tmpName);
            elseif isCOLD(tmpName)
                [sz,hash] = file.getCOLDproperties(tmpName);
            elseif isSQUID(tmpName)
                 [sz,hash] = file.getSQUIDproperties(tmpName);
            end
        end
        
        function [] = copyFile(this,target)
            sourceFile = [this.filePath this.fileName];
            copyfile(sourceFile,target);
        end
        
        % hash the name
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [hashValue] = hashName(obj)
            fName = [obj.filePath obj.fileName];
            %{
            tic
            cmd = ['echo -n ''' fName ''' | sha256sum'];
            [~,hashValue] = system(cmd);
            hashValue1 = hashValue(1:end-4);
            toc
            %}
            hashValue = hash(fName,'SHA-256','ascii');
        end
        
        % hash the name-hash and data-hash
        function [hashValue] = hashAll(obj)
            data = [obj.dataHash obj.nameHash];
            %{
            cmd = ['echo -n ''' data ''' | sha256sum'];
            [~,hashValue] = system(cmd);
            hashValue = hashValue(1:end-4);
            %}
            hashValue = hash(data,'SHA-256','ascii');
        end
        
        function [location] = storeLocation(this)
            location = getLocation(this.filePath);
        end
        
        function [CMDstring] = generateCurlCmd(this)
            CMDstring = 'curl -H "Pragma:" --speed-limit <speedLimit> --speed-time <speedTime> --retry <retry> --retry-delay <retryDelay> -o <outFile> <url>';
           
            options.speedLimit = '1024';
            options.speedTime = '30';
            options.retry = '3';
            options.retryDelay = '6';
            options.outFile = this.fileName;
           
            
            
            location = getLocation(this.filePath);
            switch location
                case 'squid'
                    urlPath = strrep(this.filePath,'squid','SQUID');
                    options.url = ['''http://proxy.chtc.wisc.edu' urlPath this.fileName ''''];
                case 'chtc'
                    fullFile = [this.filePath this.fileName];
                    options.url = ['''' file.shareS3file(fullFile) ''''];            
            end
            
            flds = fields(options);
            for e = 1:numel(flds)
                CMDstring = strrep(CMDstring,['<' flds{e} '>'],options.(flds{e}));
            end
        end
        
        function [bool] = isImage(this)
            extList = {'.tif','.jpg','.png','.bmp'};
            for e = 1:numel(this)
                [pth,nm,ext] = fileparts(this(e).fileName);
                bool(e)= any(strcmp(extList,ext));
            end
        end
        
        
    end
    
    methods (Static)
        
        % check if file is remote
        function [bool] = isRemote(fileName)
            bool = true;
            if ~isWID(fileName) && ~isIRODS(fileName)
                bool = false;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a hash link for the data-bits
        % this would allow a quick search based on the hash
        % of the file name to determine if it exists
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = makeDataLink(linkLocation,fileObject,fileName)
            % use this hash for link
            hashToUse = fileObject.dataHash;
            % use the first two for extra folder
            linkLocation = [linkLocation hashToUse(1:2) filesep];
            mmkdir(linkLocation);
            % make link for the data
            linkName = [linkLocation hashToUse];
            cmd = ['ln -s ' fileName ' ' linkName];
            system(cmd);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a hash link for the file name
        % this would allow a quick search based on the hash
        % of the file name to determine if it exists
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = makeNameLink(linkLocation,fileObject,fileName)
            hashToUse = fileObject.nameHash;
            
            linkLocation = [linkLocation hashToUse(1:2) filesep];
            mmkdir(linkLocation);
            linkName = [linkLocation hashToUse];
            cmd = ['ln -s ' fileName ' ' linkName];
            system(cmd);
        end
        
        % get size and hash values for file @ IRODS
        function [sz,hash] = getIRODSproperties(fileName)
            % generate command
            cmd = ['ils -L ' fileName];
            [~,result] = system(cmd);
            % trim result string
            result = strtrim(result);
            result = [char(10) result char(10)];
            fidx = strfind(result,char(10));

            linesPerRecord = 2;
            for e = 1:((numel(fidx)-1)/2)
                % get the first line from the eth record
                lineN = (e-1)*linesPerRecord + 1;
                line1 = result((fidx(lineN)+1):(fidx(lineN+1)-1));
                line1 = strtrim(line1);
                
                % get the second line from the eth record
                lineN = lineN+1;
                line2 = result((fidx(lineN)+1):(fidx(lineN+1)-1));
                line2 = strtrim(line2);
                
                % generate filter
                filter1 = zeros(size(line1));
                fidx1 = strfind(line1,' ');
                filter1(fidx1) = 1;
                R1 = regionprops(~logical(filter1),'PixelIdxList');
                sz{e} = line1(R1(4).PixelIdxList);
                
                % generate filter
                filter2 = zeros(size(line2));
                fidx2 = strfind(line2,' ');
                filter2(fidx2) = 1;
                R2 = regionprops(~logical(filter2),'PixelIdxList');
                hash{e} = line2(R2(1).PixelIdxList);
            end
            sz = sz{1};
            hash = hash{1};
        end
        
        % get size and hash value for file @ local
        function [sz,hash] = getLocalproperties(fileName,host)
            if nargin == 1;host = '';end
            % get size
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % list the file
            cmd = ['du -b ''' fileName ''' | cut -f1'];
            if ~isempty(host);cmd = [host '''' strrep(cmd,'''','"') ''''];end
            [out,result] = system(cmd);
            result = strtrim(result);
            sz = str2num(result);
            
            % get hash
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the hash
            cmd = ['sha256sum ''' fileName ''' | cut -f1 -d" "'];
            if ~isempty(host);cmd = [host '''' strrep(cmd,'''','"') ''''];end
            
            [~,r] = system(cmd);
            r = strtrim(r);
            hash = r;
            %fidx = strfind(r,' ');
            %hash = r(1:(fidx(1)+1));
            
        end
        
        function [sz,hash] = getCOLDproperties(fileName)
            host = ['ssh ndmiller@mir-submit.discovery.wisc.edu '];
            [sz,hash] = file.getLocalproperties(fileName,host);
        end
        
        function [sz,hash] = getSQUIDproperties(fileName)
            host = ['ssh ndmiller@submit2.chtc.wisc.edu '];
            [sz,hash] = file.getLocalproperties(fileName,host);
        end
        
        function [sz,hash] = getCHTCproperties(fileName)
            %{
            cmd = ['mc --json stat ' fileName(2:end)];
            [~,result] = system(cmd);
            result = jsondecode(result);
            %}
            cmd = ['mc --json ls ' fileName(2:end)];
            [~,result] = system(cmd);
            result = jsondecode(result);
            sz = result.size;
            hash = result.etag;
        end
        
        function [location] = shareS3file(file,expire)
            if nargin == 1;expire = [num2str(3*24) 'h '];end
            type = 'download';
            
            cmd = ['mc --json share  ' type ' --recursive -E ' expire file(2:end)];
            [~,result] = system(cmd);
            result = jsondecode(result);
            location = result.share;
        end
    end
    
end

%{
a = file('/iplant/home/nmiller/backup-Jul-11-18.tar.gz');
[sz,digest] = file.getIRODSproperties('/iplant/home/nmiller/backup-Jul-11-18.tar.gz');
b = file('/chtc/upload/test.tif');
[sz,digest] = file.getWIDproperties('/chtc/upload/test.tif');
c = file('~/what.tif');
[sz,digest] = file.getWIDproperties('/chtc/upload/test.tif');
%}
