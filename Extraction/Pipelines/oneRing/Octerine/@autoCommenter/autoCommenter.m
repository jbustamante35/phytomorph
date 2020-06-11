classdef autoCommenter < handle
    
    
    properties
        % make mirror dimension
        mirrorDimension = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/oneRing/Octerine/mirrorDimension/';
        % version of the commenter
        version = 'octerine';
        lineTAG;
        
        varBlockTag;
        descriptionBlockTag;
        acommentBlockTag;
        commentBlockStart = '%{';
        commentBlockStop = '%}';
        
        
        fileID;
    end
    
    methods
        
        
        function [obj] = autoCommenter()
            
            % make line ID tag
            obj.lineTAG.start = '<lineID>';
            obj.lineTAG.stop = autoCommenter.makeEndTag(obj.lineTAG.start);
            
            % add a version attribute to the line ID tag
            kvAttributes.version = obj.version;
            obj.lineTAG = autoCommenter.addAttribute(obj.lineTAG,kvAttributes);
            % make the line tag a comment
            obj.lineTAG.start = [' % ' obj.lineTAG.start];
            
            % make the variable block tag and its end tag
            obj.varBlockTag.start = '<varBlock>';
            obj.varBlockTag.stop = autoCommenter.makeEndTag(obj.varBlockTag.start);
        
            % make the description block
            obj.descriptionBlockTag.start = '<desBlock>';
            obj.descriptionBlockTag.stop = autoCommenter.makeEndTag(obj.varBlockTag.start);
        
            % make the com block
            obj.acommentBlockTag.start = '<comBlock>';
            obj.acommentBlockTag.stop = autoCommenter.makeEndTag(obj.acommentBlockTag.start);
        
            
        end
        
        function [b] = hasComBlock(obj,functionName)
            b = obj.functionHasXblock(functionName,obj.acommentBlockTag);
        end
        
        function [b] = hasVarBlock(obj,functionName)
            b = obj.functionHasXblock(functionName,obj.varBlockTag);
        end
        
        function [b] = hasDesBlock(obj,functionName)
            b = obj.functionHasXblock(functionName,obj.descriptionBlockTag);
        end
        
        function [] = createVarBlock(obj,functionName)
            obj.writeXblock(functionName,obj.varBlockTag);
        end
        
        function [] = createDesBlock(obj,functionName)
            obj.writeXblock(functionName,obj.descriptionBlockTag);
        end
        
        function [] = createComBlock(obj,functionName)
            aList.version = obj.version;
            tag = autoCommenter.addAttribute(obj.acommentBlockTag,aList);
            obj.writeXblock(functionName,tag);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write X block
        % write a block of commented and tagged text
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = writeXblock(obj,functionName,blockTag)
            % check if the block is present
            has_varBlock = obj.functionHasXblock(functionName,blockTag);
            % if it does not have the block
            if ~has_varBlock
                % open the file for edits
                obj.openFileForEdits(functionName);
                % comment and tag start
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',obj.commentBlockStart);
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',blockTag.start);
                fprintf(obj.fileID,'\n');
                
                % comment and tag end
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',blockTag.stop);
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',obj.commentBlockStop);
                fprintf(obj.fileID,'\n');
                % close
                obj.closeCurrentFile();
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % deterine if function has variable block
        % will scan the file and look for the block tag
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [varBlock] = functionHasXblock(obj,functionName,blockTag)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find function and open file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.openFileForRead(functionName);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % check for variable block
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            textLine = fgetl(obj.fileID);
            varBlock = false;
            varBlockLine = blockTag.start;
            while (ischar(textLine) && ~varBlock)
                varBlock = contains(textLine,varBlockLine);
                textLine = fgetl(obj.fileID);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % close file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            closeCurrentFile(obj);
        end
        
        
        function [funcPath] = openFileForEdits(obj,functionName)
            % given the function name as a string
            % find the function file
            funcPath = which(functionName);
            if isempty(funcPath)
                % assume that the functionName is the fileName
                funcPath = functionName;
            end
            % open the function for edits
            obj.fileID = fopen(funcPath,'at');
        end
        
        function [funcPath] = openFileForRead(obj,functionName)
            % find the function
            funcPath = which(functionName);
            if isempty(funcPath)
                % assume that the functionName is the fileNae
                funcPath = functionName;
            end
            obj.fileID = fopen(funcPath,'r');
        end
        
        function [] = closeCurrentFile(obj)
            fclose(obj.fileID);
        end

        function [newFile,newFunction] = interLeave(obj,functionName,weaveLine)
            
            %%%%%%%%%%%%%%%%%%%%%%
            % get set of maps
            [vec] = autoCommenter.vectorize(functionName);
            lineMap = autoCommenter.lineMap(vec);
            commentMap = autoCommenter.commentMap(vec);
            ctrlMaps = autoCommenter.ctrlMaps(vec);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%
            % open a function for read - given the functin name
            funcPath = obj.openFileForRead(functionName);
            % get the file parts
            [pth,nm,ext] = fileparts(funcPath);
            % make new file name with _interLeave as post name
            newFile = [pth filesep nm '_interLeave' ext];
            newFunction = [nm '_interLeave'];
            
            %%%%%%%%%%%%%%%%%%%%%%
            % find the line of function with the last end
            % that is used and therefore NOT in comment
            N = autoCommenter.lastEnd(funcPath);
            % open the new file for writing
            fileID_new = fopen(newFile,'w');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read and copy upto the varBlock
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            textLine = fgets(obj.fileID);
            
            %{
            writeFlag = true;
            copyFlag = true;
            inCommentFlag = false;
            isWeave = false;
            isPrior = false;
            %}
            
            lineCNT = 1;
            readCNT = 1;
            document = [];
            while (ischar(textLine))
                
                lineVec = autoCommenter.getLine(lineCNT,lineMap,vec);
                
                
                curBool = autoCommenter.isLineInMap(lineCNT,lineMap,ctrlMaps(1:2,:));
                nxtBool = autoCommenter.isLineInMap(lineCNT+1,lineMap,ctrlMaps(1:2,:));
                
                
                cmtBool = autoCommenter.isLineInMap(lineCNT,lineMap,commentMap(1,:));
                
                
                document = [document lineVec];
                
                
                
                
                %char(lineVec)
                %curBool
                %nxtBool
                
                %x = input('');
                
                if (curBool & nxtBool & ~cmtBool)
                    
                    
                    
                    
                    
                    weaveLineN = autoCommenter.lambdaCMD(weaveLine,{num2str(lineCNT)},1);
                
                    %lineID = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
                    lineID = DataHash([char(lineVec)]);
                
                    weaveLineN = [weaveLineN obj.lineTAG.start ...
                              lineID '_' num2str(lineCNT) ...
                              obj.lineTAG.stop char(10)];
                   
                          
                    
                
                    document = [document double(weaveLineN)];
                    
                end
               
                
                textLine = fgets(obj.fileID);
                readCNT = readCNT + 1;
                lineCNT = lineCNT + 1;    
                
            end
            
            fprintf(fileID_new,'%s',document);
            
            fclose(fileID_new);
            
             %   here =1 ; 
                
                %{
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % check to turn off writeFlag due to prior comments
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if (obj.isWeaveLine(textLine) && lineCNT == 1 && ~isPrior)
                    isPrior = true;
                    [textLine] = obj.stripCommentBlock(textLine);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % check to turn off writeFlag
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if autoCommenter.isStartBlockComment(textLine)
                    inCommentFlag = true;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % check to turn on writeFlag
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if (~writeFlag && autoCommenter.isStopBlockComment(textLine))
                    inCommentFlag = false;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % turn off write 
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if obj.isWeaveLine(textLine)
                    isWeave = true;
                    copyFlag = false;
                else
                    copyFlag = true;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % copy the text data into the new file
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if copyFlag
                    if ~isPrior && lineCNT == 1
                        textLine = [textLine(1:end-1) obj.lineTAG.start 'true' obj.lineTAG.stop textLine(end)];
                    end
                    fprintf(fileID_new,'%s',textLine);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % turn off write to stagger weave lines
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if (autoCommenter.isBlank(textLine) || inCommentFlag || isWeave || ...
                    autoCommenter.isComment(textLine) || isPrior || readCNT == N || ...
                    autoCommenter.isContinuedLine(textLine))
                    writeFlag = false;
                else
                    writeFlag = true;
                end
               
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % if write flag is on
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                if writeFlag
                    % 
                    lineID = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
                    lineID = DataHash([lineID lineID]);
                    
                    
                    weaveLineN = autoCommenter.lambdaCMD(weaveLine,{num2str(lineCNT)},1);
                    weaveLineN = [weaveLineN obj.lineTAG.start ...
                                    lineID '_' num2str(lineCNT) ...
                                    obj.lineTAG.stop];
                    
                    
                    fprintf(fileID_new,'%s',weaveLineN);
                    fprintf(fileID_new,'\n');
                    
                    lineCNT = lineCNT + 1;
                end
                
                
                
                
                
                
                
                %}
                % get next line
               
                
                
           % end
            
        end
        
        function [b] = isWeaveLine(obj,text)
            b = contains(text,obj.lineTAG.start);
        end
        
        function [newFuncHandle,newFunctionName,newFunctionFile] = moniterVars(obj,functionName)
            clear global
            global varLog
            varLog = containers.Map();
            
            % make mirror dimension for function
            [newFunctionName,newFunctionFile,newfList] = ...
                obj.makeMirrorDimension(functionName);
            % build the code
            cmd = autoCommenter.buildCMD('varLogger',3);
            cmd = autoCommenter.lambdaCMD(cmd,{'whos','#var1#','''f1'''});
            
            %[newFile,newFunction] = autoCommenter.copyFunctionFile(functionName,'_varLog');
            [newFunctionFile,newFunctionName] = obj.interLeave(newFunctionFile,cmd);
            [newFuncHandle] = autoCommenter.fHandleFromFile(newFunctionFile);
            
            
        end
        
        function [text] = stripCommentBlock(obj,text)
            fidx_START = strfind(text,obj.lineTAG.start);
            fidx_STOP = strfind(text,obj.lineTAG.stop);
            preText = text(1:(fidx_START(1)-1));
            postText = text(fidx_STOP(1)+numel(obj.lineTAG.stop):end);
            text = [preText postText];
        end
        
        function [newFunctionName,newFunctionFile,newfList] = makeMirrorDimension(obj,functionName)
            
            %%%%%%%%%%%%%%%%%%%%%
            mirrorID = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
            mirrorID = DataHash([mirrorID mirrorID]);
            mirrorID = mirrorID(1:10);
            % make mirror dimension co_ordinates
            mirrorPath = [obj.mirrorDimension mirrorID filesep];
            mmkdir(mirrorPath);
            
            
            %%%%%%%%%%%%%%%%%%%%%
            % add this mirror dim to the search path
            addpath(mirrorPath);
            % find default function
            funcPath = which(functionName);
            % find list of required functions
            [fList,pList] = matlab.codetools.requiredFilesAndProducts(funcPath);
            % remove recursive calls and add to list
            if ~any(strcmp(fList,funcPath))
                fList{end+1} = funcPath;
            end
            
            %%%%%%%%%%%%%%%%%%%%%
            % remove any files related to imread
            for e = 1:numel(fList)
                rm(e) = contains(fList{e},'imread');
            end
            fList(rm) = [];
            
            %%%%%%%%%%%%%%%%%%%%%
            % find the main in the list and move it to the end
            fidx = find(strcmp(fList,funcPath));
            tmp = fList{fidx};
            fList(fidx) = [];
            fList{end+1} = tmp;
            
            %%%%%%%%%%%%%%%%%%%%%
            % for each dependancy - copy into the mirror dimension
            % append the mirror dimension ID and the version of the
            % autoCommenter
            for e = 1:numel(fList)
                sourceFile = fList{e};
                [pth,nm,ext] = fileparts(sourceFile);
                
                mirrorFile = [mirrorPath nm '_' obj.version '_' mirrorID ext];
                [new_pth,new_nm,new_ext] = fileparts(mirrorFile);
                    
                oldFunctionList{e} = nm;
                newFunctionList{e} = new_nm;
                newfList{e} = mirrorFile;
                
                fprintf([sourceFile '-->' mirrorFile '\n'])
                copyfile(sourceFile,mirrorFile);
            end
            
            %%%%%%%%%%%%%%%%%%%%%
            % replace the use of the function calls to self-reference the 
            % mirror dimension
            for e = 1:numel(newfList)
                autoCommenter.replaceFunctionUse(newfList{e},oldFunctionList,newFunctionList);
            end
            
            %%%%%%%%%%%%%%%%%%%%%
            % the last in the list is the function itself
            newFunctionFile = newfList{end};
            [pth,newFunctionName,ext] = fileparts(newFunctionFile);
            
        end
        
        
       
            
 end
    
    
    methods (Static)
        
        % remove all mirron dimensions from path list
        % aka: hide the mirror dimensions
        function [] = removePath(moniker)
            if nargin == 0;moniker = 'mirrorDimension';end
            rmList = {};
            kpList = {};
            pathList = [':' path ':'];
            fidx = strfind(pathList,':');
            for e = 1:(numel(fidx)-1)
                tmp = pathList((fidx(e)+1):(fidx(e+1)-1));
                if contains(tmp,moniker)
                    rmList{end+1} = [':' tmp ':'];
                else
                    kpList{end+1} = tmp;
                end
            end
            path(kpList);
        end
        
        % given a start tag - create the proper end tag
        function [endTag] = makeEndTag(startTag)
            endTag = strrep(startTag,'<','</');
        end
        
        % copy function file with the appendValue to the file name
        function [newFile,newFunction] = copyFunctionFile(functionName,appendValue)
            % find the name of the function
            funcPath = which(functionName);
            % get the file parts
            [pth,nm,ext] = fileparts(funcPath);
            % add appendValue text to the file name
            newFile = [pth filesep nm appendValue ext];
            % make the new file name
            newFunction = [nm appendValue];
            % copy the file
            copyfile(funcPath,newFile);
        end
        
        % add attribute to a tag
        function [tag] = addAttribute(tag,aList)
            f = fields(aList);
            for e = 1:numel(f)
                tag.start = strrep(tag.start,'>',[' ' f{e} '=' aList.(f{e}) '>']);
            end
            
        end
        
        % is text the start of a block comment
        function [b] = isStartBlockComment(text)
            text = strtrim(text);
            if numel(text) >= 2
                b = strcmp(text(1:2),'%{');
            else
                b = false;
            end
        end
        
        % is text the end of a block comment
        function [b] = isStopBlockComment(text)
            text = strtrim(text);
            if numel(text) >= 2
                b = strcmp(text(1:2),'%}');
            else
                b = false;
            end
        end
        
        % is text a blank line
        function [b] = isBlank(text)
            text = strtrim(text);
            b = isempty(text);
        end
        
        % is text a commented line
        function [b] = isComment(text)
            text = strtrim(text);
            if numel(text) >= 1
                b = strcmp(text(1),'%');
            else
                b = false;
            end
        end
        
        % replace the #varN# with values
        function [cmd] = lambdaCMD(cmd,values,IDX)
            if nargin == 2;IDX = 1:numel(values);end
            if ~iscell(values);values = {values};end
            for e = 1:numel(IDX)
                cmd = strrep(cmd,['#var' num2str(IDX(e)) '#'],values{e});
            end
            
        end
        
        % input (testFunction,4) -> testFunction(#var1#,#var2#,#var3#,#var4#)
        function [cmd] = buildCMD(cmd,n)
            cmd = [cmd '('];
            for e = 1:n
               cmd = [cmd '#var' num2str(e) '#,'];
            end
            cmd(end) = [];
            cmd = [cmd ')'];
        end
        
        % compare functions
        function [] = compareFunctions(function1,function2)
            funcPath1 = which(function1);
            funcPath2 = which(function2);
            CMD = ['cmp ' funcPath1 ' ' funcPath2];
            [r,o] = system(CMD);
            if isempty(o)
                fprintf([function1 '==' function2 '\n']);
            else
                fprintf([function1 '!=' function2 '\n']);
            end
        end
        
        % ???
        function [func] = fHandleFromFile(file)
            curDir = pwd;
            [pth,nm,ext] = fileparts(file);
            cd(pth);
            func = str2func(nm);
            cd(curDir);
        end
        
        % given function file
        % replace old FunctionList with new FunctionList
        function [functionFile] = replaceFunctionUse(functionFile,oFuncList,nFuncList)
            try
                % get the file data from the original function
                [pth,nm,ext] = fileparts(functionFile);
                % create a temp function
                newFile = [pth filesep nm '_tmp' ext];
                % open the original for read
                fileID = fopen(functionFile,'r');
                % open the new for write
                newfileID = fopen(newFile,'w');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % check for variable block
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                textLine = fgets(fileID);

                while ischar(textLine)
                    % replace each in list
                    for e = 1:numel(oFuncList)
                        textLine = strrep(textLine,oFuncList{e},nFuncList{e});
                    end

                    fprintf(newfileID,'%s',textLine);
                    textLine = fgets(fileID);
                end

                fclose(fileID);
                fclose(newfileID);

                copyfile(newFile,functionFile,'f');
            catch ME
                getReport(ME);
            end
            
        end
        
        % get the number of lines in a file
        function [N] = lineCount(functionFile)
            CMD = ['cat ' functionFile '| wc -l '];
            [o,r] = system(CMD);
            r = strtrim(r);
            N = str2num(r);
        end
        
        % find the line number of the last end
        % this may not work if the end is
        % not the only text on the line
        function [N] = lastEnd(functionFile)
            % find the last end of the function
            fileID = fopen(functionFile,'r');
            % get the first line
            textLine = fgets(fileID);
            readCNT = 1;
            N = 1;
            inComment = false;
            % while can read--do
            while ischar(textLine)
                % if line is a start of block comment
                % then flag inComment = true
                if autoCommenter.isStartBlockComment(textLine)
                    inComment = true;
                end
                % if line is the end of block comment
                % then flag inComment = false
                if autoCommenter.isStopBlockComment(textLine)
                    inComment = false;
                end
                % if not in comment and line == end
                % store the last found end
                if ~inComment
                    if strcmp(strtrim(textLine),'end')
                        N = readCNT;
                    end
                end
                % get the next line
                textLine = fgets(fileID);
                % increment the line number
                readCNT = readCNT + 1;
            end
        end
        
        function [bool] = isFunctionLine(textLine)
            bool = contains(textLine,'function [');
        end
        
        function [N,lineN] = functionCount(functionFile)
            % find the last end of the function
            fileID = fopen(functionFile,'r');
            % get the first line
            textLine = fgets(fileID);
            readCNT = 1;
            lineN = 1;
            inComment = false;
            % while can read--do
            while ischar(textLine)
                % if line is a start of block comment
                % then flag inComment = true
                if autoCommenter.isStartBlockComment(textLine)
                    inComment = true;
                end
                % if line is the end of block comment
                % then flag inComment = false
                if autoCommenter.isStopBlockComment(textLine)
                    inComment = false;
                end
                % if not in comment and line == end
                % store the last found end
                if ~inComment
                    if autoCommenter.isFunctionLine(textLine)
                        lineN = [lineN readCNT];
                    end
                end
                % get the next line
                textLine = fgets(fileID);
                % increment the line number
                readCNT = readCNT + 1;
            end
            % remove the first line - dont know why
            lineN(1) = [];
            % count lineN
            N = numel(lineN);
        end
        
        function [b] = isContinuedLine(text)
            text = strtrim(text);
            if numel(text) >= 3
                b = strcmp(text(end-2:end),'...');
            else
                b = false;
            end
        end
        
        function [target] = addImageToTestProcedure(source,procedure)
            [pth,nm,ext] = fileparts(source);
            sourceDirectory = [autoCommenter.testProcedureLocation filesep procedure filesep];
            CMD = ['mkdir -p ' sourceDirectory];
            system(CMD);
            
            target =  [sourceDirectory nm ext];
            
            fprintf([source '-->' target '\n']);
            if ~exist(target,'file')
                copyfile(source,target);
            end
        end
        
        function [] = buildDepGraph(cPath)
            % start with code base as default
            if nargin == 0;cPath = '/mnt/scratch1/phytomorph_dev/';end
            % look for all the m files
            FilePath = cPath;
            FileList = {};
            FileExt = {'m'};
            FileList = fdig(FilePath,FileList,FileExt,1);

            nodes = struct([]);
            parfor e = 1:numel(FileList)
                try
                    tar = {};
                    links = struct([]);
                    [~,src] = fileparts(FileList{e});
                    nodes(e).id = src;
                    [fList,pList] = matlab.codetools.requiredFilesAndProducts(FileList{e},'toponly');
                    for k = 1:numel(fList)
                        [~,tar{k}] = fileparts(fList{k});
                    end
                    S = repmat({src},size(tar));
                    for k = 1:numel(fList)
                        links(k).source = S{k};
                        links(k).target = tar{k};
                        links(k).value = 1;
                    end
                    L{e} = links;
                    e
                catch

                end
            end
            links = struct([]);
            for e = 1:numel(L)
                links = cat(2,links,L{e});
            end


            toAdd = setdiff({links.source links.target},{nodes.id});

            for e = 1:numel(toAdd)
                more(e).id = toAdd{e};
            end
            nodes = cat(2,nodes,more);
            graphFile.nodes = nodes;
            graphFile.links = links;
            jsonStr = jsonencode(graphFile);
            jsonFileName = '~/codeGraph.json';
            fileID = fopen(jsonFileName,'w');
            fprintf(fileID,strrep(jsonStr,'\/','\\/'));
            fclose(fileID);
            

        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vectorize the function from text to double
        function [vec] = vectorize(functionFile)
            if ~isa(functionFile,'file');functionFile = which(functionFile);end
            % find the last end of the function
            fileID = fopen(functionFile,'r');
            % get the first line
            textLine = fgets(fileID);
            vec = [];
            % while can read--do
            while ischar(textLine)
                vec = [vec,double(textLine)];
                % get the next line
                textLine = fgets(fileID);
            end
            fclose(fileID);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % search the vectorized texxt for patten
        function [fvec] = find(vec,pattern)
            if ~iscell(pattern);pattern = {pattern};end
            func = @(x,y)autoCommenter.txtCmp(x,y);
            for e = 1:numel(pattern) tmpTxt = im2colF(vec,[1 numel(pattern{e})]);
                tmpTxt = autoCommenter.txtCmp(tmpTxt,double(pattern{e})');
                pad = zeros(1,numel(vec) - numel(tmpTxt));
                fvec(e,:) = [tmpTxt,pad];
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compare two vectors
        function [bool] = txtCmp(x,y)
            %if numel(y) == 1;x = x';end
            bool = all(bsxfun(@eq,x,y),1);
        end
        
        function [groupPair,groupLabel] = pairMatch(groupEnds,sz,shareEnd,vec)
            %
            groupLabel = zeros(1,size(groupEnds,2));
            
            fidxStr = find(groupEnds(1,:));
            fidxStp = find(groupEnds(2,:));
            
            % match from the last start block
            fidxStr = flip(fidxStr,2);
            for f = 1:numel(fidxStr)
                delta = fidxStp - fidxStr(f);
                delta(delta < 0) = Inf;
                [v,midx] = min(delta);
                groupPair(f,:) = [fidxStr(f),fidxStp(midx)+sz(2)-1];
                if ~shareEnd
                    fidxStp(midx) = [];
                end
                groupLabel(groupPair(f,1):groupPair(f,2)) = 1;
            end
            
            
            
            if shareEnd
                ngroupPair = [];
                len = diff(groupPair,2);
                UQE = unique(groupPair(:,2));
                for u = 1:numel(UQE)
                    fidx = find(groupPair(:,2) == UQE(u));
                    flen = len(fidx);
                    [~,midx] = max(flen);
                    ngroupPair = [ngroupPair;groupPair(fidx(midx),:)];
                end
                groupPair = ngroupPair;
            end
        end
        
        
        function [lineMap] = lineMap(vec)
            lineMap = zeros(size(vec));
            lidx = find(vec == double(char(10)));
            if lidx(end) ~= numel(vec);lidx = [lidx numel(vec)];end
            str = 1;
            for e = 1:numel(lidx)
                stp = lidx(e);
                lineMap(str:stp) = e;
                str = stp + 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make map to determine if in comment
        function [commentMap] = commentMap(vec)
            disp = false;
            % true in comment-false not
            commentMap = zeros(size(vec));
            % comment start and stop pair
            strGroupStr = {{'%{'},{'%'}};
            stpGroupStr = {{'%}'},{char(10)}};
            % flags to deterine if starts can share a stop
            shareEnd = [false,true];
            
            STR = [];
            STP = [];
            for e = 1:numel(strGroupStr)
                strChunk = [];
                for q = 1:numel(strGroupStr{e})
                    strChunk(q,:) = autoCommenter.find(vec,strGroupStr{e}{q});
                end
                strChunk = any(strChunk,1);

                stpChunk = [];
                for q = 1:numel(stpGroupStr{e})
                    stpChunk(q,:) = autoCommenter.find(vec,stpGroupStr{e}{q});
                end
                stpChunk = any(stpChunk,1);

                STR = [STR;strChunk];
                STP = [STP;stpChunk];
            end
            STR(2,:) = STR(2,:) & (~STR(1,:) & ~STP(1,:));
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % single line comments can share an end
            for e = 1:size(STR,1)
                % get the 
                strChunk = STR(e,:);
                stpChunk = STP(e,:);

                %%%%%%%%%%%%%%%%%%%%
                groupEnds = [strChunk;stpChunk];
                sz = [numel(strGroupStr{e}{1}),numel(stpGroupStr{e}{1})];
                
                %%%%%%%%%%%%%%%%%%%%
                % match blocks
                [groupPair{e},groupLabel(e,:)] = autoCommenter.pairMatch(groupEnds,sz,shareEnd(e));
                
                %%%%%%%%%%%%%%%%%%%%
                % color each pair
                tmpMap = zeros(1,numel(vec));
                for f = 1:size(groupPair{e},1)
                    if (vec(groupPair{e}(f,2)) == double(char(10)))
                        groupPair{e}(f,2) = groupPair{e}(f,2) - 1;
                    end
                    tmpMap(groupPair{e}(f,1):groupPair{e}(f,2)) = 1;
                end
                commentMap(e,:) = tmpMap;
                
                
                
                
                if disp
                    for f = 1:size(groupPair{e},1)

                        snip = char(vec(groupPair{e}(f,1):groupPair{e}(f,2)));

                        try
                            context = char(vec((groupPair{e}(f,1)-50):(groupPair{e}(f,2))+50));
                        catch ME
                            context = 'NaN';
                        end


                        fprintf(['********************************\n']);
                        %if snip(1) == double('%');snip = ['%' snip];end
                        if snip(end) == double(char(10));snip(end) = [];end
                        if context(1) == double('%');context = ['%' context];end
                        fprintf(['snip@:' strrep(snip,'%','%%') '\n']);
                        fprintf(['context@:' strrep(context,'%','%%') '\n']);
                        fprintf(['********************************\n']);

                    end
                end

            end
            commentMap = cat(1,commentMap,any(commentMap,1));
        end
        
        function [ctrlMaps] = ctrlMaps(vec,commentMap)
            % make the comment map if it does not exist
            if nargin == 1;commentMap = autoCommenter.commentMap(vec);end
            % start and end text blocks - may need to add to this list
            strGroupStr = {{'function ','classdef ','properties ','methods ','for ','if ','while ','try ','switch ',['try' char(10)]}};
            stpGroupStr = {{[char(32) 'end'],[char(9) 'end'],[char(10) 'end'],';end'}};
            % calc the length of the stop block texts - not needed
            for e = 1:numel(stpGroupStr{1})
                stpGroupLen(e) = numel(stpGroupStr{1}{e});
            end
            % loop ove the start types for the control structures
            % this should only ever be one for now
            for e = 1:numel(strGroupStr)

                % find the start blocks
                strChunk = [];
                for q = 1:numel(strGroupStr{e})
                    strChunk(q,:) = autoCommenter.find(vec,strGroupStr{e}{q});
                end
                % store each type for later filling
                strLayers = strChunk;
                strChunk = any(strChunk,1);

                % find the stop blocks
                stpChunk = [];
                for q = 1:numel(stpGroupStr{e})
                    stpChunk(q,:) = autoCommenter.find(vec,stpGroupStr{e}{q});
                end
                stpChunk = any(stpChunk,1);

                % do not look in comment blocks for code
                strChunk = strChunk & ~commentMap(3,:);
                stpChunk = stpChunk & ~commentMap(3,:);
                
                % stack for pairing process
                groupEnds = [strChunk;stpChunk];
                sz = [1 4];
                
                
                % match
                [groupPair{e},groupLabel] = autoCommenter.pairMatch(groupEnds,sz,false,vec);

                % fill the type maps
                for t = 1:size(strLayers,1)
                    locs = [ones(sum(strLayers(t,:)),1) find(strLayers(t,:))'];
                    ctrlMaps(t,:) = imfill(~logical(groupLabel),locs) & groupLabel;
                end
                ctrlMaps = cat(1,ctrlMaps,groupLabel);
                
                
                disp = 0;
                if disp
                    for f = 1:size(groupPair{e},1)
                        char(vec(groupPair{e}(f,1):groupPair{e}(f,2)))
                        %waitforbuttonpress
                    end
                end

                %{
                groupG = regionprops(logical(groupLabel(e,:)),'PixelIdxList');
                for g = 1:numel(groupG)
                    char(vec(groupG(g).PixelIdxList))
                end
                %}
            end
            %ctrlMaps = bsxfun(@times,ctrlMaps,~commentMap);
        end
        
        function [blockMaps] = blocksMaps(vec)
            % comment start and stop pair
            strGroupStr = {{'['},{'('},{'{'}};
            stpGroupStr = {{']'},{')'},{'}'}};
            % flags to deterine if starts can share a stop
            shareEnd = [false,false,false];
            % loop over search characters
            STR = [];
            STP = [];
            for e = 1:numel(strGroupStr)
                strChunk = [];
                for q = 1:numel(strGroupStr{e})
                    strChunk(q,:) = autoCommenter.find(vec,strGroupStr{e}{q});
                end
                strChunk = any(strChunk,1);

                stpChunk = [];
                for q = 1:numel(stpGroupStr{e})
                    stpChunk(q,:) = autoCommenter.find(vec,stpGroupStr{e}{q});
                end
                stpChunk = any(stpChunk,1);

                STR = [STR;strChunk];
                STP = [STP;stpChunk];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % single line comments can share an end
            for e = 1:size(STR,1)
                %%%%%%%%%%%%%%%%%%%%
                % get the or maps from above
                strChunk = STR(e,:);
                stpChunk = STP(e,:);

                %%%%%%%%%%%%%%%%%%%%
                groupEnds = [strChunk;stpChunk];
                % size feature of the start and end key words
                sz = [numel(strGroupStr{e}{1}),numel(stpGroupStr{e}{1})];
                
                %%%%%%%%%%%%%%%%%%%%
                % match blocks
                [groupPair{e},groupLabel(e,:)] = autoCommenter.pairMatch(groupEnds,sz,shareEnd(e));
                
                %%%%%%%%%%%%%%%%%%%%
                % color each pair
                tmpMap = zeros(1,numel(vec));
                for f = 1:size(groupPair{e},1)
                    % used to dis-clude the hard return 
                    % should not be needed here
                    if (vec(groupPair{e}(f,2)) == double(char(10)))
                        groupPair{e}(f,2) = groupPair{e}(f,2) - 1;
                    end
                    tmpMap(groupPair{e}(f,1):groupPair{e}(f,2)) = 1;
                end
                blockMaps(e,:) = tmpMap;
                blockLabels(e,:) = bwlabeln(tmpMap);
            end
        end
        
        function [bool] = isLineInMap(lineN,lineMap,testMap)
            for map = 1:size(testMap,1)
                bool(map) = all(testMap(map,lineMap == lineN));
            end
            bool = any(bool);
        end
        
        function [line] = getLine(lineN,lineMap,vec)
            line = vec(lineMap == lineN);
        end
        
        function [maps] = getMapSet(functionName)
            %%%%%%%%%%%%%%%%%%%%%%
            % get set of maps
            vec = autoCommenter.vectorize(functionName);
            lineMap = autoCommenter.lineMap(vec);
            commentMap = autoCommenter.commentMap(vec);
            ctrlMap = autoCommenter.ctrlMaps(vec);
            blockMaps = autoCommenter.blocksMaps(vec);
            maps = cat(1,vec,lineMap,commentMap,ctrlMap,blockMaps);
        end
        
        function [] = printMap(maps)
            for e = 1:size(maps,2)
                if maps(5,e)
                    cl = 'green';
                end
                if ~maps(5,e);cl = 'blue';end
                toPrint = strrep(char(maps(1,e)),'%','%%');
                cprintf(cl,toPrint);
            end
        end
        
        function [N,regionMap] = numRegions(map)
            [regionMap,N] = bwlabeln(map);
        end
        
        function [maps] = removeChar(maps,v)
            if ~isnumeric(v);v = double(v);end
            rm = maps(1,:) == v;
            maps(:,rm) = [];
        end
        
        function [] = stringMap(vec)
            
        end
        
        function [ORresult,ORlayer] = matchOR(vec,ORgroup)
            for q = 1:numel(ORgroup)
                ORlayer(q,:) = autoCommenter.find(vec,ORgroup{q});
            end
            % store each type for later filling
            ORresult = any(ORlayer,1);
        end
    end
    
   
end


%{


    maps = autoCommenter.getMapSet('singleCobImage_octerine');
    maps = autoCommenter.removeChar(maps,' ');
    lineN = autoCommenter.getLine(3,maps(2,:),maps(1,:));


    [blockN,blockR] = autoCommenter.numRegions(maps(end,:));

    lineR = maps(2,:);
    lineN = numel(unique(lineR));

    for e = 1:lineN
        lineV = char(maps(1,(lineR==e)));
        lineLabel = maps(end,lineR==e);
        commentValue = maps(5,lineR==e);
        blockV1 = maps(end,lineR==e);
        blockV2 = maps(end-1,lineR==e);
        blockV3 = maps(end-2,lineR==e);

        if ~any(commentValue)
            fprintf(['*******************************************\n']);
            fprintf(['line:' num2str(e) '\n']);
            fprintf(['*******************************************\n']);
            fprintf([strrep(lineV,'\n','')])
            fprintf(['-------------------------------------------\n']);
            toPrint = char(double(lineV).*blockV1);
            fprintf([strrep(toPrint,'\n','') '\n'])
            fprintf(['-------------------------------------------\n']);
            toPrint = char(double(lineV).*blockV2);
            fprintf([strrep(toPrint,'\n','') '\n'])
            fprintf(['-------------------------------------------\n']);
            toPrint = char(double(lineV).*blockV3);
            fprintf([strrep(toPrint,'\n','') '\n'])
            fprintf(['*******************************************\n']); 
        end
    end



    autoCommenter.printMap(maps)




    testFunction = 'singleCobImage_octerine';
    [newFile,newFunction] = autoCommenter.copyFunctionFile(testFunction,'_TEST');
    commenter = autoCommenter();
    commenter.makeMirrorDimension(testFunction);


    commenter.hasVarBlock(newFile)
    commenter.createComBlock(newFile)
    commenter.hasVarBlock(newFile)
    commenter.hasDesBlock(newFile)

    varFunction = commenter.moniterVars(testFunction);
    varFunction2 = commenter.moniterVars(varFunction);



    autoCommenter.compareFunctions(testFunction,varFunction2)

    [WnewFile,WnewFunction] = commenter.interLeave(newFunction,'helloWorld');

%}