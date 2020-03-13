function [List,FileList,hashList] = findObjects(searchType,parameters,toRead)
    if nargin == 2;toRead = false;end
    switch searchType
        case 'type'
            tic
            objectLocation=getenv('PHYTO_OBJECTS');
            objectLocation = [objectLocation filesep parameters];
            FilePath = objectLocation;
            FileList = {};
            FileExt = {};
            tic
            FileList = fdig(FilePath,FileList,FileExt,1);
            fprintf(['Object Search of type:' parameters ' done:' num2str(toc) '\n']);
            if nargout == 3
                for e = 1:numel(FileList)
                    [~,hashList{e}] = fileparts(FileList{e});
                end
            end
            List = FileList;
        case 'file'
            if nargin == 2;FileList={};end
            obj = fileread(parameters);
            List = jsondecode(obj);
        case 'ptrList'
            if nargin == 2;FileList={};end
            objectLocation=getenv('PHYTO_OBJECTS');
            for e = 1:numel(parameters)
                type = parameters(e).type;
                [hash] = id2hash(jsonencode(parameters(e)));
                filename = [objectLocation filesep type filesep hash(1:2) filesep hash];
                List(e) = jsondecode(fileread(filename));
                FileList{e} = filename;
            end
    end
    
    if toRead
        tic
        for e = 1:numel(FileList)
            tm = clock;
            obj = fileread(FileList{e});
            objectList(e) = jsondecode(obj);
            etm(e) = etime(clock,tm);
        end
        avergeReadTime = mean(etm);
        List = objectList;
        fprintf(['Object load of type:' parameters ':' num2str(numel(FileList)) ':' num2str(toc) ...
            '@' num2str(avergeReadTime) '\n'])
    end
    
end



function [hash] = id2hash(id)
    cmd = ['echo -n ''' id ''' | sha256sum'];
    [~,hash] = system(cmd);
    hash = hash(1:end-4);
end