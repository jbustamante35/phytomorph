classdef htcMatlabJob < htcMainShellJob
    
    
    methods
        
        function [this] = htcMatlabJob(toInit)
            this@htcMainShellJob();
            if nargin == 0;toInit = true;end
            if toInit
                % get the needed mcr
                mcrFile = htcResourceStore.getMCR();
                % ad mcr to resource list
                this.resourceList.add(mcrFile,'mcr')
            end
           
        end
        
    end
    
    
    

    
       
       
    
    
end