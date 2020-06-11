classdef htcResourceList < namedCollection
    
    properties
        
        resourceList;
        resourceName;
        resourceType;
        
        resourceActions;
    end
    
    methods
        
        function [this] = htcResourceList()
            %this@namedCollection('resourceList','cell');
            this@namedCollection('resourceList','htcResource');
        end
        
        function [] = add(this,varargin)
            % make a file a resource else assign the resource
            if isa(varargin{1},'file')
                % make file a resource
                resource = htcResource(varargin{:});
            elseif isa(varargin{1},'htcResource')
                resource = varargin{1};
            end
            % append to the end of the list
            this.data = [this.data ;resource];
            % attach the default actions based on the file type
            if isa(resource,'file.zip')
                this.resourceActions{end+1} = htcResourceList.associatedCMD(resource);
            else
                this.resourceActions{end+1} = '';
            end
        end
        
        function [] = modData(this,idx,newResource)
            if ischar(idx);idx = this.findByName(idx);end
            this.data(idx) = newResource;
        end
        
        function [] = modName(this,idx,newName)
            if ischar(idx);idx = this.findByName(idx);end
            this.data(idx) = newName;
        end
        
        function [] = remove(this,idx)
            if ischar(idx);idx = this.findByName(idx);end
            this.data(idx) = [];
        end
        
        function [idx] = findByName(this,name)
            idx = find(strcmp({this.data.name},name));
        end
        
        
        %{
        % recursive?
        function [ret] = generateCurlCmd(this,idx)
            if nargin > 1
                if ischar(idx);idx = this.findByName(idx);end
                ret = this.resourceList{idx}.generateCurlCmd();
            else
                for e = 1:numel(this.resourceList)
                    ret{e} = this.resourceList{e}.generateCurlCmd();
                end
            end
        end
        %}
    end
    
end




