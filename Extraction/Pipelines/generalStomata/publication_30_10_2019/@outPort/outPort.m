classdef outPort < handle

    properties
        
        lPath;
        rPath;
        
        filePostFix;

        fileBuffer;
    end
    
    
    methods (Abstract)
        initRemotePath(obj);
    end
    
    methods
        
        function [obj] = outPort(lPath)
            if nargin == 0; lPath = './dataPool/';postFix = '';end
            postFix = '';
            fidx = strfind(lPath,'<>');
            if ~isempty(fidx)
                postFix = lPath((fidx(1)+2):end);
                lPath = lPath(1:(fidx(1)-1));
            end
            obj.lPath = lPath;
            obj.filePostFix = postFix;
            mkdir(obj.lPath);
        end
        
        function [fileName] = localSpool(obj,fileName,data)
            %timingBlock('start','Writing data to local file.');
            [pth,nm,ext] = fileparts(fileName);
            fileName = [obj.lPath fileName];
            switch ext
                case '.csv'
                    timingBlock('note','Writing local csv file.');
                    csvwrite(fileName,data);
                case {'.tif','.jpg'}
                    timingBlock('note','Writing local image file.');
                    data = uint8(data*255);
                    imwrite(data,fileName);
                case '.mat'
                    
            end
            %timingBlock('stop');
        end
        
        function [] = remoteSpool(obj,fileName,data)
            
        end
        
        function [] = bufferSpool(obj,fileName,data)
            [fileName] = obj.localSpool(fileName,data);
            obj.fileBuffer{end+1} = fileName;
        end
        
        function [] = flushBuffer(obj,fileList)
            location = obj.rPath;
            if nargin==1;fileList = obj.fileBuffer;end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % push to irods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(location)
                % clip off ticket
                [pth,ticket,ext] = fileparts(location);
                fidx = strfind(location,'#');
                pth = location(1:(fidx(end-1)-1));
                ticket = location((fidx(end-1)+1):(fidx(end)-1));
                fprintf(['**************************************************************************\n'])
                fprintf(['Start push to irods:' num2str(numel(fileList)) ' files.\n']);tic
                fprintf(['**************************************************************************\n'])
                for e = 1:numel(fileList)
                    fprintf(['*************************************\n'])
                    fprintf(['Start push to irods:' num2str(e) ':' num2str(numel(fileList)) ' files.\n']);tic
                    fprintf(['*************************************\n'])
                    [~,tn,te] = fileparts(fileList{e});
                    targetFile = [pth filesep tn te];
                    %{
                    targetHTTP = [targetFile '#' ticket '#'];
                    targetHTTP = xform2URL({targetHTTP});
                    targetHTTP = targetHTTP{1};
                    %cmd = ['curl -X POST ' targetHTTP ' -F uploadFile=@' fileList{e}];
                    %}
                    cmd = ['iput -f -t ' ticket ' "' fileList{e} '" "' targetFile '"'];
                    fprintf(['Push Command is:' cmd '\n']);
                    [r,o] = system(cmd,'-echo');
                    %fprintf(['\npushing to file to irods:' tn '\n']);
                    fprintf(['*************************************\n'])
                    fprintf(['End push to irods:' num2str(e) ':' num2str(numel(fileList)) ' files.\n']);tic
                    fprintf(['*************************************\n'])
                end
                fprintf(['**************************************************************************\n'])
                fprintf(['End push to irods:' num2str(toc) ' seconds.\n']);
                fprintf(['**************************************************************************\n'])
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % push to irods
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % issue netwotk ticket
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fileName,iticket] = issueNetworkTicket()
            % network layer location
            fileName = ['/iplant/home/nmiller/phytoLayers'];
            % issue ticket for fileName/directory
            type = 'write';
            cmd = ['iticket create ' type ' ' '"' fileName '"'];
            [o,r] = system(cmd);
            fidx = strfind(r,':');
            iticket = r(fidx+1:end-1);
            % set expire date
            expireDate = datestr(addtodate(now,11,'day'),'YYYY-mm-dd.hh:MM:ss');
            cmd = ['iticket mod ' iticket ' expire ' expireDate];
            [o,r] = system(cmd);
            fileName = [fileName '#' iticket '#'];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % issue netwotk ticket
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fileName,iticket] = issueDataPoolTicket()
            % network layer location
            fileName = ['/iplant/home/nmiller/dataPool'];
            % issue ticket for fileName/directory
            type = 'write';
            cmd = ['iticket create ' type ' ' '"' fileName '"'];
            [o,r] = system(cmd);
            fidx = strfind(r,':');
            iticket = r(fidx+1:end-1);
            % set expire date
            expireDate = datestr(addtodate(now,11,'day'),'YYYY-mm-dd.hh:MM:ss');
            cmd = ['iticket mod ' iticket ' expire ' expireDate];
            [o,r] = system(cmd);
            fileName = [fileName '#' iticket '#'];
        end
        
        function [] = generateNewTickets()
            outPort.generateNewNetworkTicket();
            outPort.generateNewPoolTicket();
        end
        
        function [] = generateNewNetworkTicket()
            
            [rPath,iticket] = outPort.issueNetworkTicket();
            fileName = [tempdir 'networkKey.csv'];
            
            fid = fopen(fileName,'wt');
            fprintf(fid, rPath);
            fclose(fid);
            
            cmd = ['iput -f ' fileName ' /iplant/home/nmiller/publicData/networkKey.csv'];
            [r,o] = system(cmd);
        end
        
        function [] = generateNewPoolTicket()
            
            [rPath,iticket] = outPort.issueDataPoolTicket();
            fileName = [tempdir 'poolKey.csv'];
            
            fid = fopen(fileName,'wt');
            fprintf(fid, rPath);
            fclose(fid);
            
            cmd = ['iput -f ' fileName ' /iplant/home/nmiller/publicData/poolKey.csv'];
            [r,o] = system(cmd);
        end
        
        function [fileName] = formatFileName(nameStruct,delimiter)
            if nargin == 1;delimiter='_';end
            fExt = nameStruct.ext;
            f = fields(nameStruct);
            f = setdiff(f,'ext');
            fileName = '';
            for e = 1:numel(f)
                fileName = [fileName '{' f{e} delimiter nameStruct.(f{e}) '}'];
            end
            fileName = [fileName '.' fExt];
        end

    end


    

end