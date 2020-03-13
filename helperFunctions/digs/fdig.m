function [FileList] = fdig(FilePath,FileList,FileExt,verbose)
try
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
    
    [r,o] = system(CMD);
    oidx = strfind(o,char(10));
    oidx = [0 oidx];
    for e = 1:(numel(oidx)-1)
        FileList{e} = o(oidx(e)+1:oidx(e+1)-1);
    end
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% directory of input FilePath
    dirList = dir(FilePath);
    ridx = strcmp({dirList.name},'.') | strcmp({dirList.name},'..');
    dirList(ridx) = [];
    %%% directory of input FilePath
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if not empty
    if size(dirList,2) ~= 0
        % for each object in the list
        for listing = 1:size(dirList,1)
            %%% init the next path
            current_Path = [FilePath dirList(listing).name];            
            typed_path = regexprep(current_Path,filesep,[filesep filesep]);
            typed_path = regexprep(typed_path,'%',['%%']);
            % if directory
            if dirList(listing).isdir
                % report 
                if verbose
                    fprintf(['Looking at:' typed_path '\n']);
                end
                % call gdig
                FileList = gdig([current_Path filesep],FileList,FileExt,verbose);
            % if file
            else
                % look for the .
                tidx = strfind(dirList(listing).name,'.');
                % if not empty and a type in the list - add
                if ~isempty(tidx)                    
                    if any(strcmp(FileExt,dirList(listing).name(tidx(end)+1:end))) || any(strcmp(FileExt,'*'))
                        FileList{end+1} = current_Path;
                    end
                end
            end
        end
    end
    %}
catch ME
    ME
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
