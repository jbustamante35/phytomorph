classdef htcArgumentStore < objectFSstore
    
    
    properties (Constant)
        defaultDataPath = '/mnt/scratch1/htcArguments/';
        defaultDatabaseFile = '/mnt/scratch1/arguments.db'
    end
    
    
    properties
        database;
    end
    
    methods
        
        function [this] = htcArgumentStore(basePath)
            if nargin == 0;basePath = htcArgumentStore.defaultDataPath;end
            this@objectFSstore(basePath);
            this.database = htcArgumentStore.defaultDatabaseFile;
            mksqlite('open',this.database);
            mksqlite(['CREATE TABLE IF NOT EXISTS arguments '...
                '(id INTEGER PRIMARY KEY AUTOINCREMENT, type TEXT,hash TEXT,uuid TEXT)']);
            mksqlite(['CREATE TABLE IF NOT EXISTS functions '...
                '(id INTEGER PRIMARY KEY AUTOINCREMENT, functionName TEXT,uuid TEXT)']);
        end
        
        
        function [fileName] = add(this,argument)
            % insert into database
            mksqlite('open',this.database);
            % store the argument type not the container type
            type = class(argument.value);
            % capture the hash
            hash = argument.hash;
            % capture the uuid
            uuid = argument.uuid;
            % build the insert command
            sql = ['INSERT INTO arguments VALUES (?,?,?,?)'];
            % insert the arguements
            mksqlite(sql,{'',type,hash,uuid});
            % close the database
            mksqlite('close');
            % write to disk
            fileName = this.write(argument,'mat');
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
        
    end
    
    
end