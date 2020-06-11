classdef htcResourceStore < objectFSstore
    
    properties (Constant)
        defaultPath = '/mnt/scratch1/htcResources/';
    end
    
    properties
        database;
    end
    
    methods
        
        function [this] = htcResourceStore(basePath)
            if nargin == 0;basePath = htcResourceStore.defaultPath;end
            this@objectFSstore(basePath);
            this.database = [basePath 'resource.db'];
            mksqlite('open',this.database);
            mksqlite(['CREATE TABLE IF NOT EXISTS resources '...
                '(id INTEGER PRIMARY KEY AUTOINCREMENT,name TEXT, type TEXT,hash TEXT)']);
        end
        
        function [] = add(this,resource,name)
            % insert into database
            mksqlite('open',this.database);
            type = class(resource);
            hash = resource.uuid;
            sql = ['INSERT INTO resources VALUES (?,?,?,?)'];
            mksqlite(sql,{'',name,type,hash});
            mksqlite('close');
            
            % write to disk
            this.write(resource);
        end
        
        function [resource] = find(this,queryValue,queryType)
            mksqlite('open',this.database);
            if nargin < 3;queryType = 'name';end
            if strcmp(queryValue(1),'?')
                queryValue = strrep(queryValue,'?',queryType);
                sql = ['SELECT name,type,hash FROM resources WHERE REGEX(' queryValue ') NOT NULL'];
                results = mksqlite(sql);
            else
                sql = ['SELECT name,type,hash FROM resources WHERE ' queryType '=''' queryValue ''''];
                results = mksqlite(sql);
            end
            mksqlite('close');
            % find the file
            resource = find@objectFSstore(this,'file',results.hash);    
        end
        
        function [resourceTable] = list(this)
            mksqlite('open',this.database);
            sql = ['SELECT * FROM resources'];
            results = mksqlite(sql);
            resourceTable = table;
            for e = 1:numel(results)
                resourceTable.name{e} = results(e).name;
                resourceTable.type{e} = results(e).type;
            end
            mksqlite('close');
            resourceTable
        end
        
      
    end
    
    methods (Static)
         function [fl] = getMCR(versionSTR)
             % connect to the default store
            store = htcResourceStore();
            if nargin == 0
                versionSTR = version;
                fidx = strfind(versionSTR,'.');
                versionSTR = versionSTR(1:fidx(3));
                versionSTR = strrep(versionSTR,'.','');
            end
            sql = ['?,"\[mcr=' versionSTR '\]"'];
            fl = store.find(sql);
        end
    end
end