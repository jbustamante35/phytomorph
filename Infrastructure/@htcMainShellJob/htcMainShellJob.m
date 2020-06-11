classdef htcMainShellJob < shellScript & htcJob
    
    
    properties (Constant)
        % report the os type and bit level
        mainline01 = 'uname -m';
        mainline02 = 'uname';
        
        % init the irods environment
        mainline001 = 'export IRODS_ENVIRONMENT_FILE=$PWD/irods_environment.json';
        mainline002 = 'export ENVIRONMENT_VAR_HOME=$PWD';
        mainline003 = 'export IRODS_AUTHENTICATION_FILE=$PWD/pwfile';
    end
    
    methods
        
        
        function [this] = htcMainShellJob()
            this@htcJob;
            this@shellScript;
        end
        
        
    end
    
    
    
    
    

    
       
       
    
    
end