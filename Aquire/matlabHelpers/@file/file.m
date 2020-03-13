classdef file < doid

    properties
        localPath;
        fileName;
        size;
        
        dataHash;
        nameHash;
        fullHash;
        
        stage = 1;
        
    end

    methods
        function [obj] = file(fName)
            obj = obj@doid();
            tmp = dir(fName);
            obj.size = tmp.bytes;
            [obj.localPath,nm,ext] = fileparts(fName);
            obj.localPath = [obj.localPath filesep];
            obj.fileName = [nm ext];
            obj.dataHash = obj.hashData();
            obj.nameHash = obj.hashName();
            obj.fullHash = obj.hashAll();
            % overwrite the uuid as the full-hash=hash(hash(data)+hash(name))
            obj.uuid = obj.fullHash;
        end
        
        
        function [hash] = hashData(obj)
            cmd = ['sha256sum ' obj.localPath obj.fileName];
            [~,r] = system(cmd);
            fidx = strfind(r,' ');
            hash = r(1:(fidx(1)+1));
        end
        
        function [hash] = hashName(obj)
            fName = [obj.localPath obj.fileName];
            cmd = ['echo -n ''' fName ''' | sha256sum'];
            [~,hash] = system(cmd);
            hash = hash(1:end-4);
        end
        
        function [hash] = hashAll(obj)
            data = [obj.dataHash obj.nameHash];
            cmd = ['echo -n ''' data ''' | sha256sum'];
            [~,hash] = system(cmd);
            hash = hash(1:end-4);
        end
        
    end
    
    methods (Static)
        
        function [] = makeDataLink(linkLocation,fileObject,fileName)
            hashToUse = fileObject.dataHash;
            
            linkLocation = [linkLocation hashToUse(1:2) filesep];
            mmkdir(linkLocation);
            linkName = [linkLocation hashToUse];
            cmd = ['ln -s ' fileName ' ' linkName];
            system(cmd);
        end
        
        function [] = makeNameLink(linkLocation,fileObject,fileName)
            hashToUse = fileObject.nameHash;
            
            linkLocation = [linkLocation hashToUse(1:2) filesep];
            mmkdir(linkLocation);
            linkName = [linkLocation hashToUse];
            cmd = ['ln -s ' fileName ' ' linkName];
            system(cmd);
        end
        
        
        
    end
    
end
