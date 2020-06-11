classdef htcResource < namedCollection
    
    properties
        rType;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [this] = htcResource(value,name,type)
            % this is the typle of name,value,type
            this@namedCollection(name,value);
            % if the type is not given then infer the type
            if nargin < 3
                % infer type
                type = htcResource.inferType(value);
            end
            % set the tpe
            this.rType = type;
        end
    end
    
    methods (Static)
        % infer the type based on this static functio
        function [cmd] = associatedCMD(file,type)
            if nargin == 1;type = inferType(file);end

            cmd = '';
            [~,~,ext] = fileparts(file.fileName);

            switch type
                case 'file.zip'
                    cmd = shellCMD('unzip ');
                case 'file.exe'
                    cmd = shellCMD('chmod +x ');
                case 'file.gz'
                    cmd = shellCMD('tar xvf ');
            end
        end
        
        
        function [type] = inferType(resource)
            if isa(resource,'file')
                [~,~,ext] = fileparts(resource);
                if ~isempty(ext)
                    type = [class(resource) ext];
                else
                    fprintf(['Unable to infer type for file.\n']);
                    fprintf(['file:' file.filePath file.fileName]);
                end
            elseif isa(resource,'htcResourceList')
                 type = class(resource);
            end
        end
        
    end
end