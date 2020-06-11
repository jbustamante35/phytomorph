function [FileList,metaData] = fdig(FilePath,FileList,FileExt,verbose,metaDataFlag)
    try
        if nargin < 5;metaDataFlag = false;end
        metaData = struct('dataHash','','size','');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prepare the file extension
        if ~isempty(FileExt)
            CMD = ['find ''' FilePath ''' -type f \( '];
            for e = 1:numel(FileExt)
               CMD = [CMD '-name \*.' FileExt{e} ' -o '];
            end
            CMD(end-2:end) = [];
            CMD = [CMD '\)'];
        else
             CMD = ['find ''' FilePath ''' -type f'];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isCOLD(FilePath)
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % gather the basics
            CMD = strrep(CMD,'''','"');
            baseCMD = CMD;
            % add the checksum if needed
            if metaDataFlag
                CMD = [baseCMD ' | xargs dataHash '];
            end
            CMD = ['ssh ndmiller@submit2.chtc.wisc.edu ''' CMD ''''];
            [r,o] = system(CMD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parse
            oidx = strfind(o,char(10));
            oidx = [0 oidx];
            for e = 1:(numel(oidx)-1)
                line = o(oidx(e)+1:oidx(e+1)-1);
                
                if ~isempty(line)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if meta data is present
                    if metaDataFlag
                        sidx = strfind(line,' ');
                        metaData(numel(FileList)+1).dataHash = line(1:(sidx(1)-1));
                        line = line((sidx(1)+2):end);
                    end
                    FileList{end+1} = line;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if metaData requested then grab the size
            if metaDataFlag
                % add the checksum if needed
                CMD = [baseCMD ' | xargs du -b '];
                CMD = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' CMD ''''];
                [r,o] = system(CMD);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % parse
                oidx = strfind(o,char(10));
                oidx = [0 oidx];
                for e = 1:(numel(oidx)-1)
                    line = o(oidx(e)+1:oidx(e+1)-1);
                    sidx = strfind(line,char(9));
                    byteValue = line(1:(sidx(1)-1));
                    fileName = line((sidx(1)+1):end);
                    idx = find(strcmp(FileList,fileName));
                    metaData(idx).size = byteValue;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    CMD = strrep(CMD,'''','"');
    baseCMD = CMD;
    if metaData
        CMD = [baseCMD ' | xargs dataHash '];
    end
    CMD = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' CMD ''''];
    [r,o] = system(CMD);
    %}
            
        elseif isSQUID(FilePath)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % gather the basics
            CMD = strrep(CMD,'''','"');
            baseCMD = CMD;
            % add the checksum if needed
            if metaDataFlag
                CMD = [baseCMD ' | xargs sha256sum '];
            end
            CMD = ['ssh ndmiller@submit2.chtc.wisc.edu ''' CMD ''''];
            [r,o] = system(CMD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parse
            oidx = strfind(o,char(10));
            oidx = [0 oidx];
            for e = 1:(numel(oidx)-1)
                line = o(oidx(e)+1:oidx(e+1)-1);
                if ~isempty(line)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if meta data is present
                    if metaDataFlag
                        sidx = strfind(line,' ');
                        metaData(numel(FileList)+1).dataHash = line(1:(sidx(1)-1));
                        line = line((sidx(1)+2):end);
                    end
                    FileList{end+1} = line;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if metaData requested then grab the size
            if metaDataFlag
                % add the checksum if needed
                CMD = [baseCMD ' | xargs du -b '];
                CMD = ['ssh ndmiller@submit2.chtc.wisc.edu ''' CMD ''''];
                [r,o] = system(CMD);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % parse
                oidx = strfind(o,char(10));
                oidx = [0 oidx];
                for e = 1:(numel(oidx)-1)
                    line = o(oidx(e)+1:oidx(e+1)-1);
                    sidx = strfind(line,char(9));
                    byteValue = line(1:(sidx(1)-1));
                    fileName = line((sidx(1)+1):end);
                    idx = find(strcmp(FileList,fileName));
                    metaData(idx).size = byteValue;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        elseif isLOCAL(FilePath)
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % gather the basics
            CMD = strrep(CMD,'''','"');
            baseCMD = CMD;
            % add the checksum if needed
            if metaDataFlag
                CMD = [baseCMD ' | xargs -I {} sha256sum "{}" '];
            end
            %CMD = ['ssh ndmiller@submit2.chtc.wisc.edu ''' CMD ''''];
            [r,o] = system(CMD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parse
            oidx = strfind(o,char(10));
            oidx = [0 oidx];
            for e = 1:(numel(oidx)-1)
                line = o(oidx(e)+1:oidx(e+1)-1);
                if ~isempty(line)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if meta data is present
                    if metaDataFlag
                        sidx = strfind(line,' ');
                        metaData(numel(FileList)+1).dataHash = line(1:(sidx(1)-1));
                        line = line((sidx(1)+2):end);
                    end
                    FileList{end+1} = line;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if metaData requested then grab the size
            if metaDataFlag
                % add the checksum if needed
                CMD = [baseCMD ' | xargs -I {} du -b "{}" '];
                %CMD = ['ssh ndmiller@submit2.chtc.wisc.edu ''' CMD ''''];
                [r,o] = system(CMD);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % parse
                oidx = strfind(o,char(10));
                oidx = [0 oidx];
                for e = 1:(numel(oidx)-1)
                    line = o(oidx(e)+1:oidx(e+1)-1);
                    sidx = strfind(line,char(9));
                    byteValue = line(1:(sidx(1)-1));
                    fileName = line((sidx(1)+1):end);
                    idx = find(strcmp(FileList,fileName));
                    metaData(idx).size = byteValue;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %{
            [r,o] = system(CMD);
            oidx = strfind(o,char(10));
            oidx = [0 oidx];
            for e = 1:(numel(oidx)-1)
                line = o(oidx(e)+1:oidx(e+1)-1);
                FileList{e} = line;
            end
            %}
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
       

        
        
        

    catch ME
        ME
    end
    
end

%{
%%%
% Useful examples
%%%
% ALL DATA
FilePath = 'W:\';
FilePath = '/mnt/spaldingdata/nate/mirror_images/rue/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);

FilePath = '/home/nate/iplant/tassels/2015/';
FileList = {};
FileExt = {'jpg'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc

FilePath = '/mnt/tetra/nate/fixPOP/next/20180221_Rack2_Camera6/';
FileList = {};
FileExt = {'tiff'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc

FilePath = '/mnt/spaldingdata/nate/';
FileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FileList = fdig(FilePath,FileList,FileExt,1);
w = find(contains(FileList,'seed'))

%}
