function [FileList] = dig(FilePath,FileList,FileExt,verbose,toFoid)
    if nargin < 2;FileList = {};end
    if nargin < 3;FileExt = {};end
    if nargin < 4;verbose = false;end
    if nargin < 5;toFoid = false;end
    if toFoid;global store;end

    if nargin == 3;verbose=false;end
    location = getLocation(FilePath);
    
    if ~iscell(FileExt);FileExt={FileExt};end
    
    fprintf(['search@' location '@' FilePath '\n']);
    switch location
        case 'irods'
            [FileList] = idig(FilePath,FileList,FileExt,verbose);
        case 'cold'
            [FileList,metaData] = fdig(FilePath,FileList,FileExt,verbose,toFoid);
        case 'local'
            [FileList,metaData] = fdig(FilePath,FileList,FileExt,verbose,toFoid);
        case 'chtc'
            [FileList] = wdig(FilePath,FileList,FileExt,verbose,toFoid);
        case 'squid'
            [FileList,metaData] = fdig(FilePath,FileList,FileExt,verbose,toFoid);
    end
    
    
    if toFoid
        if isa(store,'objectFSstore')
            tmpFileList = {};
            for e = 1:numel(FileList)
                tic
                tmpFileList{e} = file(FileList{e});
                fprintf(['converted:' num2str(e) ':' num2str(numel(FileList)) ':' num2str(toc) '\n']);
            end
        else
            tmpFileList = {};
            for e = 1:numel(FileList)
                %{
                if strcmp(location,'squid')
                    N = [15 20];
                    P = randi(N(2)) + N(1);
                    fprintf(['Pausing for SQUID.\n']);
                    for p = 1:P
                        fprintf('.');
                        pause(1);
                    end
                    fprintf('\n');
                end
                %}
                
                tic
                tmpFileList{e} = file(FileList{e},metaData(e));
                fprintf(['converted:' num2str(e) ':' num2str(numel(FileList)) ':' num2str(toc) '\n']);
            end
        end
        FileList = tmpFileList;
    end
    
    
end
