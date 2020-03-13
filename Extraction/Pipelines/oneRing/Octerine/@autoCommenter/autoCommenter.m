classdef autoCommenter < handle
    
    
    properties
        
        mirrorDimension = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/oneRing/Octerine/mirrorDimension/';
        
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
            
            obj.lineTAG.start = '<lineID>';
            obj.lineTAG.stop = autoCommenter.makeEndTag(obj.lineTAG.start);
            
            
            kvAttributes.version = obj.version;
            obj.lineTAG = autoCommenter.addAttribute(obj.lineTAG,kvAttributes);
            obj.lineTAG.start = [' % ' obj.lineTAG.start];
            
            
            obj.varBlockTag.start = '<varBlock>';
            obj.varBlockTag.stop = autoCommenter.makeEndTag(obj.varBlockTag.start);
        
        
            obj.descriptionBlockTag.start = '<desBlock>';
            obj.descriptionBlockTag.stop = autoCommenter.makeEndTag(obj.varBlockTag.start);
        
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
        % write variable block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = writeXblock(obj,functionName,blockTag)
            has_varBlock = obj.functionHasXblock(functionName,blockTag);
            
            
            if ~has_varBlock
               
                obj.openFileForEdits(functionName);

                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',obj.commentBlockStart);
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',blockTag.start);
                fprintf(obj.fileID,'\n');
                
                


                
                
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',blockTag.stop);
                fprintf(obj.fileID,'\n');
                fprintf(obj.fileID,'%s',obj.commentBlockStop);
                fprintf(obj.fileID,'\n');

                obj.closeCurrentFile();
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % deterine if function has variable block
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
            funcPath = which(functionName);
            if isempty(funcPath)
                % assume that the functionName is the fileNae
                funcPath = functionName;
            end
            obj.fileID = fopen(funcPath,'at');
        end
        
        function [funcPath] = openFileForRead(obj,functionName)
            
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

        function [newFile,newFunction] = interLeave(obj,functionName,weaveLine,type)
            
            funcPath = obj.openFileForRead(functionName);
            
            
            [pth,nm,ext] = fileparts(funcPath);
            newFile = [pth filesep nm '_interLeave' ext];
            newFunction = [nm '_interLeave'];
            
            
            N = autoCommenter.lastEnd(funcPath);
            
            fileID_new = fopen(newFile,'w');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read and copy upto the varBlock
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            textLine = fgets(obj.fileID);
            writeFlag = true;
            copyFlag = true;
            inCommentFlag = false;
            isWeave = false;
            isPrior = false;
            lineCNT = 1;
            readCNT = 1;
            while (ischar(textLine))
                
             
                
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
                
                
                % get next line
                textLine = fgets(obj.fileID);
                readCNT = readCNT + 1;
                
                
            end
            
        end
        
        function [b] = isWeaveLine(obj,text)
            b = contains(text,obj.lineTAG.start);
        end
        
        function [newFuncHandle,newFunctionName,newFunctionFile] = moniterVars(obj,functionName)
            clear global
            
            global varLog
            varLog = containers.Map();
            
             [newFunctionName,newFunctionFile,newfList] = ...
                 obj.makeMirrorDimension(functionName);
            
          
            
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
            
            mirrorID = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
            mirrorID = DataHash([mirrorID mirrorID]);
            mirrorPath = [obj.mirrorDimension filesep mirrorID filesep];
            mkdir(mirrorPath);
            
            addpath(mirrorPath)
            
            funcPath = which(functionName);
            
            
            
            
            [fList,pList] = matlab.codetools.requiredFilesAndProducts(funcPath);
            
            if ~any(strcmp(fList,funcPath))
                fList{end+1} = funcPath;
            end
            
            
            
            fidx = find(strcmp(fList,funcPath));
            tmp = fList{fidx};
            fList(fidx) = [];
            fList{end+1} = tmp;
            
            
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
            
            
            
            for e = 1:numel(newfList)
                autoCommenter.replaceFunctionUse(newfList{e},oldFunctionList,newFunctionList);
            end
            
            newFunctionFile = newfList{end};
            [pth,newFunctionName,ext] = fileparts(newFunctionFile);
            
        end
        
        
        
 end
    
    
    methods (Static)
        function [endTag] = makeEndTag(startTag)
            endTag = strrep(startTag,'<','</');
        end
        
        function [newFile,newFunction] = copyFunctionFile(functionName,appendValue)
            funcPath = which(functionName);
            [pth,nm,ext] = fileparts(funcPath);
            newFile = [pth filesep nm appendValue ext];
            newFunction = [nm appendValue];
            copyfile(funcPath,newFile);
        end
        
        function [tag] = addAttribute(tag,aList)
            f = fields(aList);
            for e = 1:numel(f)
                tag.start = strrep(tag.start,'>',[' ' f{e} '=' aList.(f{e}) '>']);
            end
            
        end
        
        function [b] = isStartBlockComment(text)
            text = strtrim(text);
            if numel(text) >= 2
                b = strcmp(text(1:2),'%{');
            else
                b = false;
            end
        end
        
        function [b] = isStopBlockComment(text)
            text = strtrim(text);
            if numel(text) >= 2
                b = strcmp(text(1:2),'%}');
            else
                b = false;
            end
        end
        
        function [b] = isBlank(text)
            text = strtrim(text);
            b = isempty(text);
        end
        
        function [b] = isComment(text)
            text = strtrim(text);
            if numel(text) >= 1
                b = strcmp(text(1),'%');
            else
                b = false;
            end
        end
        
        function [cmd] = lambdaCMD(cmd,values,IDX)
            if nargin == 2;IDX = 1:numel(values);end
            if ~iscell(values);values = {values};end
            for e = 1:numel(IDX)
                cmd = strrep(cmd,['#var' num2str(IDX(e)) '#'],values{e});
            end
            
        end
        
        function [cmd] = buildCMD(cmd,n)
            cmd = [cmd '('];
            for e = 1:n
               cmd = [cmd '#var' num2str(e) '#,'];
            end
            cmd(end) = [];
            cmd = [cmd ')'];
        end
        
        
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
        
        function [func] = fHandleFromFile(file)
            curDir = pwd;
            [pth,nm,ext] = fileparts(file);
            cd(pth);
            func = str2func(nm);
            cd(curDir);
        end
        
        function [functionFile] = replaceFunctionUse(functionFile,oFuncList,nFuncList)
            [pth,nm,ext] = fileparts(functionFile);
            newFile = [pth filesep nm '_tmp' ext];
            fileID = fopen(functionFile,'r');
            newfileID = fopen(newFile,'w');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % check for variable block
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            textLine = fgets(fileID);
            
            while ischar(textLine)
                
                for e = 1:numel(oFuncList)
                    textLine = strrep(textLine,oFuncList{e},nFuncList{e});
                end
                
                fprintf(newfileID,'%s',textLine);
                textLine = fgets(fileID);
            end
            
            fclose(fileID);
            fclose(newfileID);
            copyfile(newFile,functionFile);
            
        end
        
        function [N] = lineCount(fileName)
            CMD = ['cat ' fileName '| wc -l '];
            [o,r] = system(CMD);
            r = strtrim(r);
            N = str2num(r);
        end
        
        function [N] = lastEnd(functionFile)
            fileID = fopen(functionFile,'r');
            textLine = fgets(fileID);
            readCNT = 1;
            N = 1;
            inComment = false;
            while ischar(textLine)
                
                if autoCommenter.isStartBlockComment(textLine)
                    inComment = true;
                end
                
                if autoCommenter.isStopBlockComment(textLine)
                    inComment = false;
                end
                
                if ~inComment
                    if strcmp(strtrim(textLine),'end')
                        N = readCNT;
                    end
                end
                
                textLine = fgets(fileID);
                readCNT = readCNT + 1;
            end
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
        if nargin == 0;cPath = '/mnt/scratch1/phytomorph_dev/';end
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

    end
    
   
end


%{


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