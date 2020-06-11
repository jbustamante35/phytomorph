classdef htcBag < doid
    
    properties
        % dag representation of bag
        dag;
        % list of resource files for setting up the environment to run
        resourceList;
        % mcr
        MCR_version = 'v717';
    end
    
    methods
        
        function [this] = htcBag()
            this.resourceList = htcResourceList();
            this.dag = htcDag();
        end
        
    end
    
end