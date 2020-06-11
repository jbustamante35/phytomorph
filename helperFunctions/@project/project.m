classdef project < oid
    
    
    properties (Constant)
        defaultLocation = '/mnt/scratch1/phytomorph_dev/projectStore/';
        defaultOutputLocation = '/mnt/spaldingdata/nate/octerineDataStore/results/';
    end
    
    properties
        %
        name;
        openedFiles;
        activeDirectory;
        
        dataCollections;
        functionList;
        resultsCollections = {};
    end
    
    methods
    
        function [this] = project(name,openedFiles,activeDirectory)
            if nargin < 1;name = 'default';end
            if nargin < 2;openedFiles = {};end
            if nargin < 3;activeDirectory = pwd;end
            this.name = name;
            this.openedFiles = openedFiles;
            this.activeDirectory = activeDirectory;
        end
        
        % print the project
        function [] = print(this,leaderTxt,color,verbose)
            if nargin < 3;color = 'black';end
            if nargin < 4;verbose = false;end
            cprintf(color,[leaderTxt this.name '\n']);
        end
      
        function [b] = isActive(this)
            % get the bname if this were active
            bName = this.generateBname();
            b = exist(bName,'file');
        end
        
        function [bName] = generateBname(this,location)
            if nargin == 1;location = project.defaultLocation;end
            sessionPID = num2str(feature('getpid'));
            bName = [location this.name filesep '{projectName_' this.name '}{sessionID_' sessionPID '}.bin'];
        end
        
        function [] = writeActiveFile(this)
            activeFile = project.findActive();
            if isempty(activeFile);activeFile = generateBname(this);end
            % get current/active name
            activeName = project.getActiveName(activeFile);
            if ~strcmp(activeName,this.name)
                activeFile = project.moveActiveFile(activeFile,this.name);
            end
            activeTime = datestr(now,'yyyy-MM-dd hh:mm:ss');
            fileID = fopen(activeFile,'a');
            fprintf(fileID,activeTime);
            fclose(fileID);
        end
        
        function [] = activate(this)
            this.writeActiveFile();
            this.loadFiles();
            this.activateDirectory();
        end
        
        function [] = loadFiles(this)
            % list of opened files
            filenames = this.openedFiles;
            % open the files
            for ii = 1:length(filenames)
                if exist(filenames{ii}, 'file')
                    matlab.desktop.editor.openDocument(filenames{ii});
                else
                    warning(['File "' filenames{ii} '" was not found']);
                end
            end
        end
        
        function [] = activateDirectory(this)
            try
                cd(this.activeDirectory);
            catch
                warning(['Directory "' this.activeDirectory '" does not exist'])
            end
        end
        
        function [returnLocation] = generateReturnLocation(this)
            returnDate = datestr(now,'yyyy-MM-dd hh:mm:ss');
            returnDate = strrep(returnDate,'-','_');
            returnDate = strrep(returnDate,':','_');
            returnDate = strrep(returnDate,' ','_');
            returnLocation = [project.defaultOutputLocation ...
                            this.name filesep 'returns' filesep returnDate filesep];
        end
            
        function [] = profileFunction(this,funcName,collectionName,N)
            % get requested function
            for e = 1:numel(this.functionList)
                fidx(e) = strcmp(this.functionList{e}.name,funcName);
            end
            fidx = find(fidx);
            func = this.functionList{fidx};
            
            % get requested collection
            for e = 1:numel(this.dataCollections)
                fidx(e) = strcmp(this.dataCollections{e}.name,collectionName);
            end
            fidx = find(fidx);
            collection = this.dataCollections{fidx};
            
            % run profile and save resutlts
            func.profile(collection,N);
            
            % save this
            project.save(this);
        end
        
        
            
    end
    
    methods (Static)
        
        function [prj] = import(oldType)
            prj = project(oldType.ProjectName,oldType.OpenedFiles,oldType.ActiveDir);
        end
        
        function [varargout] = load(projectName)
            % close the current project
            project.close(false);
            % get current project
            prj = project.getProject(projectName);
            % activate the projet
            prj.activate();
            % assign output
            if nargout == 1;varargout{1} = prj;end
        end
        
        function [prj] = getProject(projectName)
            % build the requested project name
            fName = [project.defaultLocation projectName filesep projectName '.mat'];
            prj = load(fName);
            prj = prj.this;
        end
        
        % get the active project
        function [prj] = getCurrent()
            [activeFile,activeName,activeMat] = project.findActive();
            prj = load(activeMat);
            prj = prj.this;
        end
        
        
        function [] = save(projectName,baseLocation)
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            if nargin < 2;baseLocation = project.defaultLocation;end
          
            
            
            
            % get opened documents
            openDocuments = matlab.desktop.editor.getAll;
            % make the list of opened files
            openedFiles = {openDocuments.Filename};
            % get active directory
            activeDirectory = pwd();
            
            % save current project
            if nargin == 0
                
                [activeFile,activeName,activeMat] = project.findActive();
                if ~isempty(activeMat)
                    prj = load(activeMat);prj = prj.this;
                    prj.activeDirectory = activeDirectory;
                    prj.openedFiles = openedFiles;
                    project.save(prj);
                else
                    fprintf(['No active project. Can not save.\n']);
                end
               
            
               
            % save new project
            elseif nargin == 1
                % save the project
                if isa(projectName,'project')
                    % make this
                    this = projectName;
                    % make save location
                    location = [baseLocation this.name filesep];
                    % make directory
                    mmkdir(location);
                    % make mar file
                    fName = [location this.name '.mat'];
                    % save
                    save(fName,'this');
                    % report
                    disp(['Project "' projectName.name '" was saved']);
                else
                    % make new project
                    prj = project(projectName,openedFiles,activeDirectory);
                    % save project
                    project.save(prj);
                    % load project
                    project.load(projectName);
                end
            end
            
           
        end
      
        % list all the projects and highlight the current
        % for the current, display the Y=F(X).
        function [] = list()
            % get the list of project files
            FilePath = project.defaultLocation;
            FileList = {};
            FileExt = {'mat'};
            FileList = fdig(FilePath,FileList,FileExt,0);

            % make the list
            fprintf('---------------------------\n');
            fprintf('List of available projects:\n')
            fprintf('---------------------------\n');
            for e = 1:numel(FileList)
                prj = load(FileList{e},'this');prj = prj.this;
                ident = '\t';color = 'black';
                
               
                
                if prj.isActive()
                    
                    
                    ident = '------->';color = '*red';
                    leaderTxt = [ident num2str(e) '.'];
                    prj.print(leaderTxt,color);
                
                  
                    
                    spacer = ['\t   |-=========================\n'];
                    cprintf('black',spacer);
                    
                    
                    data_ident = ['\t   |-'];
                    dataColor = '*blue';
                    cprintf(dataColor,[data_ident '-----------------\n']);
                    cprintf(dataColor,[data_ident 'X-Collections:\n'])
                    cprintf(dataColor,[data_ident '-----------------\n']);
                    if ~isempty(prj.dataCollections)
                        for d = 1:numel(prj.dataCollections)
                            tmpD = [data_ident num2str(d) '.' prj.dataCollections{d}.name '\n'];
                            cprintf(dataColor,tmpD);
                        end
                    else
                        for d = 1:1
                            tmpD = [data_ident '\n'];
                            cprintf(dataColor,tmpD);
                        end
                    end
                    cprintf(dataColor,[data_ident '-----------------\n']);
                    
                    
                    
                    spacer = ['\t   |-=========================\n'];
                    cprintf('black',spacer);
                    
                    
                    
                    
                    
                    % function list
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    function_ident = ['\t   |-'];
                    functionColor = '*green';
                    %functionColor = '*[0.1,.25,0.18]';
                    %functionColor = '*green';
                    cprintf(functionColor,[function_ident '-----------------\n']);
                    cprintf(functionColor,[function_ident 'F-Collections:\n'])
                    cprintf(functionColor,[function_ident '-----------------\n']);
                    if ~isempty(prj.functionList)
                        for d = 1:numel(prj.functionList)
                            functionColor = '*[0,1,0]';
                            % body indent
                            body_ident = ['\t      |-'];
                            % signature
                            signature = prj.functionList{d}.getSignature();
                            
                            % prepare main print statement header for func
                            tmpD = [function_ident num2str(d) '.' prj.functionList{d}.name '\n'];
                            %tmpD = [tmpD '<-----' profileString  '\n'];
                            cprintf(functionColor,tmpD);
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % report properties
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % get the last profile
                            [~,profileStruct] = prj.functionList{d}.lastProfile;
                            % get function name
                            body = prj.functionList{d}.getFunctionBody();
                            % add to report
                            profileStruct.body = ['Function Body: ' body];
                            % report signature
                            profileStruct.signature = ['Signature:' signature.string];
                            
                            flds = fields(profileStruct);
                            
                            fidx = find(strcmp(flds,'signature'));
                            
                            flds = {flds{fidx} flds{:}};
                            flds(fidx+1) = [];
                            
                            % add units to timing
                            profileStruct.(flds{3}) = [profileStruct.(flds{3}) ' sec'];
                            
                            for f = 1:numel(flds)
                                 tmpD = [body_ident profileStruct.(flds{f}) '\n'];
                                 cprintf(functionColor,tmpD);
                            end
                            
                            
                        end
                        
                    
                    else
                        for d = 1:1
                            tmpD = [function_ident '\n'];
                            cprintf(functionColor,tmpD);
                        end
                    end
                    cprintf(functionColor,[function_ident '-----------------\n']);
                    
                    

                    spacer = ['\t   |-=========================\n'];
                    cprintf('black',spacer);
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % results list
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    results_ident = ['\t   |-'];
                    resultsColor = '*magenta';
                    cprintf(resultsColor,[results_ident '-----------------\n']);
                    cprintf(resultsColor,[results_ident 'Y-Collections:\n'])
                    cprintf(resultsColor,[results_ident '-----------------\n']);
                    
                    
                    if ~isempty(prj.resultsCollections)
                        for d = 1:numel(prj.resultsCollections)
                            % body indent
                            body_ident = ['\t      |-'];
                            
                            % prepare main print statement header for func
                            tmpD = [results_ident num2str(d) '.' prj.resultsCollections{d}.name '\n'];
                            %tmpD = [tmpD '<-----' profileString  '\n'];
                            cprintf(resultsColor,tmpD);
                            
                            
                            
                            % report properties
                            %%%%%%%%%%%%%%%%%%%%%
                            %location = 
                            percentFail = numel(prj.resultsCollections{d}.failIndex)/...
                                    numel(prj.resultsCollections{d}.inputDataCollection.FileList);
                            numberFail = numel(prj.resultsCollections{d}.failIndex);
                            propertiesStruct.failPercent = ['Percent Fail:' num2str(percentFail)];
                            propertiesStruct.numberFail = ['Number Fail:' num2str(numberFail)];
                            flds = fields(propertiesStruct);
                            for f = 1:numel(flds)
                                 tmpD = [body_ident propertiesStruct.(flds{f}) '\n'];
                                 cprintf(resultsColor,tmpD);
                            end
                            
                            extra = 1;
                            
                            
                        end
                        
                    
                    else
                        for d = 1:1
                            tmpD = [function_ident '\n'];
                            cprintf(resultsColor,tmpD);
                        end
                    end
                    
                    
                    
                    cprintf(resultsColor,[function_ident '-----------------\n']);
                    
                    spacer = ['\t   |-=========================\n'];
                    cprintf('black',spacer);
                    
                    
                    
                else
                    leaderTxt = [ident num2str(e) '.'];
                    prj.print(leaderTxt,color);
                end
               
            end
        end
        
        function [] = close(openDefault)
            if nargin == 0;openDefault=true;end
            
            % check if current is modified
            [b,activeName] = project.isModified();
            % gather response
            if b
                ans = '';
                while isempty(ans)
                    ans = input(['save:' activeName '(y,N)?\n'],'s');
                    if isempty(ans)
                        fprintf(['respond!\n']);
                    end
                end
                
                % save if requested
                if strcmp(lower(ans),'y')
                    project.save();
                end
            end
            % close all documents
            openDocuments = matlab.desktop.editor.getAll;
            openDocuments.close;
            % open default if needed
            if openDefault;project.load('default');end
            
        end
          
        function [activeFile,activeName,activeMat] = findActive()
            activeName = '';activeMat = '';
            % get session ID
            sessionPID = num2str(feature('getpid'));
            % look for all active files
            FilePath = project.defaultLocation;
            FileList = {};
            FileExt = {'bin'};
            FileList = gdig(FilePath,FileList,FileExt,false);
            % search for active files containing session ID
            fidx = contains(FileList,sessionPID);
            activeFile = FileList(fidx);
            
            if ~isempty(activeFile)
                activeFile = activeFile{1};
                [~,nm] = fileparts(activeFile);
                [kvp] = findKVP(nm,'projectName');
                activeName = getValue(kvp);
                activeMat = [project.defaultLocation  activeName filesep activeName '.mat'];
            end
        end
        
        function [b,activeName] = isModified()
            b = false;
            [activeFile,activeName,activeMat] = project.findActive();
            if ~isempty(activeMat)
                % load the active
                this = load(activeMat);
                this = this.this;
                % get list of open documents
                openDocuments = matlab.desktop.editor.getAll;
                fn_opened = {openDocuments.Filename};
                % compare lists
                if length(this.openedFiles) ~= length(fn_opened)
                    b = true;
                else
                    for ii=1:length(fn_opened)
                        if ~strcmpi(fn_opened{ii},this.openedFiles{ii})
                            b = true;
                        end
                    end
                end
            end
           
        end
        
        
        % collection API - this is the X
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [] = attachDataCollection(name,location,varargin)
            newCollection = dataCollection(name,location,varargin{:});
            curProject = project.getCurrent();
            curProject.dataCollections{end+1} = newCollection;
            project.save(curProject);
        end
      
        function [] = removeCollection(name)
            fidx = project.searchForCollection(name);
            if ~isempty(fidx)
                prj = project.getCurrent();
                prj.dataCollections(fidx) = [];
                project.save(prj);
            end
        end
        
        function [fidx,collection] = searchForCollection(name)
            prj = project.getCurrent();
            for e = 1:numel(prj.dataCollections)
                fidx(e) = strcmp(prj.dataCollections{e}.name,name);
            end
            fidx = find(fidx);
            collection = prj.dataCollections{fidx};
        end
        
        function [] = renameCollection(name,newName)
            fidx = project.searchForCollection(name);
            if ~isempty(fidx)
                prj = project.getCurrent();
                prj.dataCollections{fidx}.name = newName;
                project.save(prj);
            end
        end
        
        function [collection] = getCollection(name)
            [~,collection] = project.searchForCollection(name);
        end
        
        function [] = refreshCollection(name)
            prj = project.getCurrent();
            
            if strcmp(name,'all')
                for e = 1:numel(prj.dataCollections)
                    prj.dataCollections{e}.refreshList();
                end
            else
                fidx = project.searchForCollection(name);
                if ~isempty(fidx)
                    prj.dataCollections{fidx}.refreshList();
                end
            end
          
            
            project.save(prj);
        end
        
        % function API - this is the F
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [] = attachFunction(name,func)
            func = formalFunc(name,func);
            curProject = project.getCurrent();
            curProject.functionList{end+1} = func;
            project.save(curProject);
        end
        
        function [func] = getFunction(name)
            [~,func] = project.searchForFunction(name);
        end
        
        function [fidx,func] = searchForFunction(name)
            prj = project.getCurrent();
            func = '';
            if ~isempty(prj.functionList)
                for e = 1:numel(prj.functionList)
                    fidx(e) = strcmp(getName(prj.functionList{e}),name);
                end
                fidx = find(fidx);
                func = prj.functionList{fidx};
            end
        end
        
        function [] = deleteFunction(name)
            if ischar(name);fidx = project.searchForFunction(name);
            else; fidx = name;end
            curProject = project.getCurrent();
            curProject.functionList(fidx) = [];
            project.save(curProject);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [sz] = collectionSize(collectionName)
            curProject = project.getCurrent();
            collection = project.getCollection(collectionName);
            sz = numel(collection.FileList);
        end
        
        % eval/verb API  - this is F(X) - perhaps the = sign?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % run a collection (X) against a function (F)
        % totalN is the number to run
        function [] = runCollection(collectionName,totalN,functionName)
            if nargin < 3;functionName = 'main';end
            
            % get the current project
            curProject = project.getCurrent();
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the X
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the named collection
            collection = project.getCollection(collectionName);
            % get N run all collection if not given N
            if nargin < 2;totalN = numel(collection.FileList);end
            if totalN == inf;totalN = numel(collection.FileList);end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the F
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the named function
            func = project.getFunction(functionName);
            % get the active project name
            projectName =  curProject.name;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make the Y
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make return collection
            returnCollection = resultsContainer(projectName,functionName);
            % get the return location for the return collection
            returnLocation = returnCollection.path;


            % loop and run the data
            strTime = clock;
            % get random set
            collection = collection.getRandomSet(totalN);
            % number of processors
            numP = 12;
            
            FileList = collection.FileList;
            
            for e = 1:totalN
                % start the clock for the eth eval of F
                i_strTime{e} = clock;
                % get the eth file
                fileName = FileList{e};
                % call the function
                results{e} = func.operate(fileName,returnLocation);
                % get the stop time
                i_stpTime{e} = clock;
                % get the stop-start for the eth job
                i_deltaTime(e) = etime(i_stpTime{e},i_strTime{e});
                % get the average job time
                meanJobTime = mean(i_deltaTime);
                % calc the expected amount of time for all jobs
                expectedTotal = totalN*meanJobTime;
                % calc the time spent on jobs
                i_wallTime = etime(i_stpTime{e},strTime);
                % calc the fraction done
                fraction = i_wallTime / expectedTotal;
                % calc the remaining time
                remainTime = expectedTotal - i_wallTime;
                %%%%%%%%%%%%%%%%%%%
                %{
                % get the total time elapsed since starting running jobs
                i_wallTime = etime(i_stpTime,strTime);
                % 
                % calc the number done 
                i_done = i_wallTime / i_deltaTime;
                % 
                i_done = i_done * numP;
                fraction = i_done / totalN;
                totalTime = (totalN*i_deltaTime)/numP;
                remainTime = (1-fraction)*totalTime;
                %}
                %%%%%%%%%%%%%%%%%%%
                cprintf('red',['~percent:' num2str(fraction) '~remain:' num2str(remainTime) '(s)\n'])
            end
            %stpTime = clock;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % look for fails
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for e = 1:numel(results)
                failIndex(e) = isa(results{e},'MException');
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % config the Y to know the fail index
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set the input collection and the fail index
            returnCollection.inputDataCollection = collection;
            returnCollection.failIndex = find(failIndex);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save results
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % store the return collection
            curProject.resultsCollections{end+1} = returnCollection;
            % save the project
            project.save(curProject);
        end
        
        function [results] = runDatum(collectionName,imageName,functionName)
            if nargin < 3;functionName = 'main';end
            func = project.getFunction(functionName);
            curProject = project.getCurrent();
            collection = project.getCollection(collectionName);
            fileList = collection.searchByName(imageName);
            returnLocation = curProject.generateReturnLocation();
            mmkdir(returnLocation);
            for e = 1:numel(fileList)
                results{e} = func.operate(fileList{e},returnLocation);
            end
        end
        
        % this is th Y API
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [] = clean(functionName)
            % trash location
            trashLocation = '/mnt/spaldingdata/nate/octerineDataStore/results/trash/';
            % default to main
            if nargin < 1;functionName = '*';end
            % get the current project
            curProject = project.getCurrent();
            % assign functionList
            if ~isa(functionName,'cell')
                if ~strcmp(functionName,'*')
                    toCleanList = {functionName};
                else
                    toCleanList = curProject.functionList;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            for i = 1:numel(toCleanList)
                functionName = toCleanList{i}.name;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % scan the default output location 
                oPath = [resultsContainer.defaultOutputLocation ...
                    '{projectName_' curProject.name '}' filesep ...
                    '{functionName_' functionName '}' filesep];
                % look in (project,function)
                cdir = dir(oPath);cdir(1:2) = [];
                % keep only the directories
                cdir(~[cdir.isdir]) = [];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % glue to make a list of runs
                for e = 1:numel(cdir)
                    runList_total{e} = [oPath cdir(e).name filesep];
                end
                % make a list of successful runs
                for e = 1:numel(curProject.resultsCollections)
                    runList_good{e} = curProject.resultsCollections{e}.path;
                end
                % remove the good to get the fail
                runList_fail = setdiff(runList_total,runList_good);
                % move each fail folder
                for e = 1:numel(runList_fail)
                    movefile(runList_fail{e},trashLocation);
                end
            end
                
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [value] = getActiveName(activeFile)
            [~,nm] = fileparts(activeFile);
            kvp = findKVP(nm,'projectName');
            value = getValue(kvp);
        end
        
        function [targetFile] = moveActiveFile(activeFile,targetProjectName)
            [~,nm] = fileparts(activeFile);
            kvp = findKVP(nm,'projectName');
            value = getValue(kvp);
            newKVP = strrep(kvp,value,targetProjectName);
            targetFile = strrep(activeFile,kvp,newKVP);
            targetFile = strrep(targetFile,value,targetProjectName);
            movefile(activeFile,targetFile);
        end
        
        % profile a function (default is main)
        % use the collection = collectionname
        % use N (default is 10)
        function [] = profile(collectionName,N,functionName)
            if nargin < 2;N = 10;end
            if nargin < 3;functionName = 'main';end
            % get the current project
            curProject = project.getCurrent();
            curProject.profileFunction(functionName,collectionName,N);
        end
        
    end
    
end






