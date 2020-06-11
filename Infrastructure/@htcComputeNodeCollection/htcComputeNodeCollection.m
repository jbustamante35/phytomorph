classdef htcComputeNodeCollection < namedCollection
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is a named collection of htcComputeNodes
    % a collection of compute nodes will share 
    % 1-resources
    % 2-hardware requirements
    % 3-submit file
    %%%%%
    % compute nodes belong to collections
    % jobs run on compute nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        
        % properties list
        props;
        % plus-properties list
        pprops;
        
        % resourceList
        resourceList;
    end
    
    methods
        %
        function [this] = htcComputeNodeCollection(name)
            % this is a named collection of compute nodes
            this@namedCollection(name,'htcComputeNode');
            
            % init the resource list for this container
            this.resourceList = htcResourceList();
            
            
            %%%%%%%%%%%%%%%%%%%%
            % init the props for the hardware and config
            % the props are held in a struct container
            %%%%%%%%%%%%%%%%%%%%
            % init the default properties
            this.props = struct();
            this.initHardwareProps();
            this.initConfigProps();
            %
            this.pprops = struct();
            this.pprops.Group = '"spalding"';
        end
        
        % default properties for a colletion of compute nodes
        function [] = initHardwareProps(this)
            this.props.universe = 'vanilla';
            this.props.request_memory = 8000;
            this.props.cpus = 1;
            this.props.disk = 8000000;
        end
        % default config properties
        function [] = initConfigProps(this)
            this.props.priority = 9;
            this.props.should_trander_files = 'YES';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [] = generateSUBMITfile(this,oPath)
            oFile = [oPath this.name '.submit'];
            fileID = fopen(oFile,'w');
            flds = fields(this.props);
            for e = 1:numel(flds)
                value = this.props.(flds{e});
                if isnumeric(value);value = num2str(value);end
                fprintf(fileID,'%s\n',[flds{e} '=' value]);
            end
            flds = fields(this.pprops);
            for e = 1:numel(flds)
                value = this.pprops.(flds{e});
                if isnumeric(value);value = num2str(value);end
                fprintf(fileID,'%s\n',[flds{e} '=' value]);
            end
            fclose(fileID);
        end
        
        
    end
end