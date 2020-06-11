classdef objectFSstore < handle
    
    properties
        % base path for the data store
        basePath;
        % base path for softlinks listed via data hash
        base_dataSoftLinks;
        % base path for softlinks listed via name hash
        base_nameSoftLinks;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % file object store constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = objectFSstore(basePath)
            % get the environmental variable for the phytoobjects
            if nargin == 0;basePath = getenv('PHYTO_OBJECTS');end
            % set the base path
            obj.basePath = basePath;
            % mkdir basepath
            mmkdir(obj.basePath);
            % set the data hash softlinks location
            obj.base_dataSoftLinks = [obj.basePath 'file' filesep 'dataHash' filesep];
            mmkdir(obj.base_dataSoftLinks);
            % set the name hash softlinks location
            obj.base_nameSoftLinks = [obj.basePath 'file' filesep 'nameHash' filesep];
            mmkdir(obj.base_nameSoftLinks);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate an object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [object] = generate(obj,type,varargin)
            % call the constructure of 'type' with varargin
            [input] = obj.buildVarInput(varargin{:});
            % build constructor
            constructor = str2func(['@' input type input]);
            % call constructor
            object = constructor(varargin{:});
            % write the object to disk
            fileName = obj.write(object);
            
            % special instructions when making a file object
            if isa(object,'file')
                % cross list the file under the data hash
                file.makeDataLink(obj.base_dataSoftLinks,object,fileName);
                % cross list the file under the name hash
                file.makeNameLink(obj.base_nameSoftLinks,object,fileName);
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % invoke function on object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = invoke(obj,object,func,varargin)
            % attach object to input
            varargin = {object,varargin{:}};
            % build input algebra
            [input] = obj.buildVarInput(varargin{:});
            % build function call
            func = str2func(['@' input func input]);
            % invoke
            func(varargin{:});
            % write the object
            obj.write(object);
            % loop over each input to see if its oid and save
            for e = 1:numel(varargin)
                if isa(varargin{e},'oid')
                    obj.write(varargin{e});
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [object] = read(obj,otr)
            % if not a char then must be a dptr or higher oid
            if ~isa(otr,'char')
                if isa(otr,'dptr');otr = otr.dereference();end
                [otr] = buildFileName(obj,otr);
            else
                [otr] = buildFileName(obj,otr);
            end
            % read the data from file
            object = fileread(otr);
            % convert from json to obect
            object = oid.fromJSON(object);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fileName] = write(obj,otw,writeType)
            % write the type to disk
            if nargin < 3;writeType = 'json';end
            % copy the object before writing
            otw = copy(otw);
            % detach the copy
            otw.detach();
            % build the file name
            [fileName] = buildFileName(obj,otw);
            [pth,nm,ext] = fileparts(fileName);
            % make output locatin
            mmkdir(pth);
            switch writeType
                case 'json'
                    % open the file to write
                    fileID = fopen(fileName,'w');
                    % encode the json object
                    % this needs to be better and tested
                    jsonotw = jsonencode(otw);
                    % write the json daata to disk
                    fprintf(fileID,'%s',jsonotw);
                    % close the file
                    fclose(fileID);
                case 'mat'
                    fileName = [fileName '.mat'];
                    save(fileName,'otw');
            end

            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find objects of type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ret] = find(obj,type,uuid)
            if nargin == 2
                FilePath = [obj.basePath filesep type];
                FileList = {};
                FileExt = {};
                tic
                FileList = fdig(FilePath,FileList,FileExt,1);
                ret = FileList;
                fprintf(['Object Search of type:' parameters ' done:' num2str(toc) '\n']);
            elseif nargin == 3
                otr = oid(type,uuid);
                ret = obj.read(otr);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reattach object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = reattach(obj,object,loadList)
            if nargin == 2;loadList.uuid = {};loadList.obj = {};end
            prop = properties(object);
            for p = 1:numel(prop)
                if isa(object.(prop{p}),'oid')
                    % loop over the oid array and convert to ptr
                    for e = 1:numel(object.(prop{p}))
                        if isa(object.(prop{p})(e),'dptr')
                            % check if object has been loaded -
                            toLoadPtr = object.(prop{p})(e);
                            toLoadUuid = toLoadPtr.refs.uuid;
                            if ~isempty(setdiff(toLoadUuid,loadList.uuid))
                                % read the object
                                object.(prop{p})(e) = obj.read(object.(prop{p})(e));
                                loadList.uuid{end+1} = toLoadUuid;
                                loadList.obj{end+1} = object.(prop{p})(e);
                                % here is recursion
                                obj.reattach(object.(prop{p})(e),loadList);
                            else
                                fidx = strcmp(loadList.uuid,toLoadUuid);
                                object.(prop{p})(e) = loadList.obj{fidx};
                            end
                            
                        end
                    end
                    
                end
            end
        end
    end
    
    methods (Access=private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is the single point for objects to have file names
        % attached to them
        % note: may 26, 2020
        % filename <- [base - type - hash[1:2] - hash]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [fileName] = buildFileName(obj,object)
            % this line was traded for the uuid - why
            % hash = object.hash();
            % this line was traded for a case statement
            % hash = object.uuid;
            % put back to uuid then overwritten in file constructor
            
            
            % mod to load from hash char
            if ~isa(object,'char')
                hash = object.uuid;
            else
                hash = object;
            end
            
            fileName = [obj.basePath object.type filesep hash(1:2) filesep hash];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build an input string for varargin
        % repeating units of x -> (x1,x2,...,xn)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [input] = buildVarInput(obj,varargin)
            N = numel(varargin);
            base = 'x';
            input = '(';
            for e = 1:N
                input = [input base num2str(e) ','];
            end
            
            if numel(input) == 1
                input = [input ')'];
            else
                input(end) = ')';
            end
        end
        
        
    end
end