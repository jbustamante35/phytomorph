classdef dataCollection < doid
    
    properties (Constant)
        ImageExt = {'tiff','tif','jpg','TIF','TIFF','PNG','JPG'};
    end
    
    properties
        name;
        FileList;
        baseLocation;
        FileExt;
        createDate;
        isSet;
        isFrozen;
        isFile;
        szBytes = NaN;
    end
    
    
    methods
        
        % name of data collection
        function [this] = dataCollection(name,baseLocation,FileExt,isSet,isFile)
            % make default to be no name
            if nargin == 0;name = '';baseLocation='';end
            % default is a file ext of images
            if nargin <= 2;FileExt = dataCollection.ImageExt;end
            %
            if nargin < 4;isSet = false;end
            if nargin < 5;isFile = false;end
            if isempty(FileExt);FileExt = dataCollection.ImageExt;end
            % if the first argument is a char - then name
            if isa(name,'char')
                % if baseLocation is char - then computable list
                if isa(baseLocation,'char')
                    this.name = name;
                    this.baseLocation = baseLocation;
                    this.FileExt = FileExt;
                    this.createDate = datestr(now,'yyyy-MM-dd hh:mm:ss');
                    this.isSet = isSet;
                    this.isFrozen = false;
                    this.isFile = isFile;
                    this.refreshList();
                % if baseLocation is cell - then frozen - non-computable
                elseif isa(baseLocation,'cell')
                    this.name = name;
                    this.FileList = baseLocation;
                    this.baseLocation = '';
                    this.FileExt = '';
                    this.createDate = datestr(now,'yyyy-MM-dd hh:mm:ss');
                    % if cell or not
                    this.isSet = isa(baseLocation{1},'cell');
                    this.isFrozen = true;
                    this.isFile = isFile;
                end
            % if first arg is a cell - then anonymous filelist
            elseif isa(name,'cell')
                this.FileList = name;
                this.name = '';
                this.baseLocation = '';
                this.FileExt = {};
                this.createDate = datestr(now,'yyyy-MM-dd hh:mm:ss');
                % if cell or not
                this.isSet = isa(name{1},'cell');
                this.isFrozen = true;
                this.isFile = isFile;
            end
        end
        
        function [] = refreshList(this)
            if ~this.isFrozen
                if ~this.isSet
                    this.FileList = dig(this.baseLocation,{},this.FileExt,1,this.isFile);
                else
                    this.FileList = sdig(this.baseLocation,{},this.FileExt,1,this.isFile);
                end
                if this.isFile
                    this.szBytes = 0;
                    for e = 1:numel(this.FileList)
                        this.szBytes = this.szBytes + this.FileList{e}.size;
                    end
                end
            end
        end
        
        function [ret] = searchByName(this,name)
            fidx = find(contains(lower(this.FileList),name));
            ret = this.FileList(fidx);
        end
        
        function [ret] = getRandomSet(this,N)
            if (N == numel(this.FileList))
                ret = this;
            else
                ridx = randi(numel(this.FileList),N,1);
                subSet = this.FileList(ridx);
                collectionName = ['randomSubSet@' this.name];
                ret = dataCollection(collectionName,subSet);
            end
          
        end
        %{
        function [n] = numel(this)
            n = numel(this.FileList);
        end
        %}
    end
    
    
    methods (Static)
        
        function [resultCollection] = buildFromDirectory(source)
            [FileList] = dig(source,{},{},false,true);
            resultCollection = dataCollection.groupByStatus(FileList);
        end
        
        function [resultCollection] = groupByStatus(FileList)
            for e = 1:numel(FileList)
                [pth,nm,ext] = fileparts(FileList{e});
                isStatus(e) = contains(nm,'{status_');
            end
            fidx = find(isStatus);
            for e = 1:numel(fidx)
                [pth,nm,ext] = fileparts(FileList{fidx(e)});
                [kvp] = findKVP(nm,'fileName');
                groupKey{e} = getValue(kvp);
                [kvp] = findKVP(nm,'status');
                statusKey{e} = getValue(kvp);
            end
            for e = 1:numel(FileList)
                [pth,nmList{e},ext] = fileparts(FileList{e});
            end
            for e = 1:numel(groupKey)
                fidx = contains(nmList,groupKey{e});
                tmpSub = FileList(fidx);
                resultCollection(e) = resultsCollection(tmpSub);
                resultCollection(e).setStatus(statusKey{e});
            end
        end
        
       
        
    end
    
end