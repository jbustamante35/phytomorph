classdef partialFunction
    properties
        func;
        algorithm;
    end
    
    methods
        function [obj] = partialFunction(func,varargin)
            if ~isa(func,'function_handle')
                tempdir = [tempname filesep];
                mkdir(tempdir);
                matFile = [tempdir func '.mat'];
                CMD = ['iget -f /iplant/home/nmiller/publicData/' func '.mat ' matFile];
                
                system(CMD);
                load(matFile);
            else
                obj.func = func;
                obj.algorithm = varargin{1};
            end
        end
        
        function [] = save(obj,filename)
            save(filename,'obj','-v7.3');
        end
        
        function [] = publish(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this functionality has been documented on 10.30.19
            % this functionality has been upgraded on 10.30.19 
            % 1) perform backup pre upload so that i can roll back to other
            % functions if i brea things
            % 2) publish on squid too (NOT DONE YET)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % render partial function to mat file under algorithm name
            tempdir = [tempname filesep];
            tempdir = ['/mnt/tetra/nate/' tempdir];
            % make local locatin for the function
            mkdir(tempdir);
            % name of the mat file 
            matFile = [tempdir obj.algorithm '.mat'];
            % save the mat file local
            obj.save(matFile);
             % get the file parts
            [pth nm ext] = fileparts(matFile);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make the backup directory
            CMD = ['imkdir -p /iplant/home/nmiller/publicData/functionBackups/' nm filesep];
            % make the backup directory
            system(CMD);
            fprintf(['Making backup location for:' nm '\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            backupName = [nm '_backup_' datestr(datetime,'ss_hh_dd_mm_yyyy') ext];
            % comamnd to copy the old function to the backup
            CMD = ['icp /iplant/home/nmiller/publicData/' nm ext ' '...
                    '/iplant/home/nmiller/publicData/functionBackups/' backupName];
            system(CMD);
            fprintf(['Backing up:' nm '-->' backupName '\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % construct command to save the mat file
            CMD = ['iput -f ' matFile ' /iplant/home/nmiller/publicData/'];
            % render the mat file to the cyverse store
            system(CMD);
            fprintf(['Deposting new copy.\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % execute system commands to make the file downloadable
            CMD = ['ichmod read anonymous /iplant/home/nmiller/publicData/' nm ext];
            system(CMD);
            CMD = ['curl -H ''Cache-Control: no-cache'' https://data.cyverse.org/dav-anon/iplant/home/nmiller/publicData/' nm ext];
            system(CMD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function [o] = evalFunction(obj,X)
            o = obj.func(X);
        end

        function [w] = editArgs(obj)
            w = functions(obj.func);
            w = w.workspace{1};
        end

        
    end
end

%{
    J = 3;
    g = @(X)X*J;
    f = partialFunction(g,'test');
    f.publish();
    h = partialFunction('test');
    l1 = f.evalFunction(4)
    l2 = h.evalFunction(4)
%}