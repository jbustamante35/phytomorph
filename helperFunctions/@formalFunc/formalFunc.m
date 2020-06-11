classdef formalFunc < doidmm
    
    properties (Constant)
        storeLocation = '/mnt/myCompile/octerineFlow/';
    end
    
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % human name of function
        name;
        % func
        func;
        % version
        version;
        % signature
        signature;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % profile table
        outputProfileTable;
        % state
        computeState;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods
        
        function [this] = formalFunc(name,func)
            this@doidmm();
            if nargin == 1;func=name;name='';end
            % init the name
            this.name = name;
            % default compute state
            this.computeState = true;
            if ~isa(func,'function_handle')
                this.func = funcFromName(func);
                %signature = formalFunc.signatureFromName(func);
                %this.func = str2func(['@' signature.pinIN signature.body signature.pinIN]);
            else
                % init the function - this is a function handle
                %this.func = func;
                % func from name
                this.func = funcFromName(functionName);
            end
            % init the version
            this.version = 0;
            % init the signature
            this.resetSignature();
            % init the profile table
            this.initProfileTable();
        end
        
        function [n] = numArgumentsFromSubscript(obj,s,ic)
            n = 0;
        end
        
        % subs ref for API (),{} calls to the base function
        function [varargout] = subsref(this,subs)
            switch subs(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',this,subs);
                case '()'
                    if this.computeState
                        in = subs.subs;
                        if (nargin == 1)
                            if isa(in{1},'argfile')
                                in = in{1}.load();
                                [varargout{1:nargout}] = this.func(in{:});
                            else
                                %%%%%%%%%%%%%%%%%%%%%%%%%
                                % normal compute
                                %%%%%%%%%%%%%%%%%%%%%%%%%
                                % get signature
                                sig = this.signature;
                                % check for in-dim lower
                                % if true - try to invoke with lower dim
                                % signature
                                if (sig.size(2) ~= numel(in))
                                    tmpIN = formalFunc.removeEndInputs(sig.pinIN,sig.size(2)-nargin);
                                    tmpFunc = str2func(['@' tmpIN sig.body tmpIN]);
                                    [varargout{1:nargout}] = tmpFunc(in{:});
                                else
                                    [varargout{1:nargout}] = this.func(in{:});
                                end
                            end
                        else
                            %%%%%%%%%%%%%%%%%%%%%%%%%
                            % normal compute
                            %%%%%%%%%%%%%%%%%%%%%%%%%
                            % get signature
                            sig = this.signature;
                            % check for in-dim lower
                            % if true - try to invoke with lower dim
                            % signature
                            if (sig.size(2) ~= numel(in))
                                tmpIN = formalFunc.removeEndInputs(sig.pinIN,sig.size(2)-nargin);
                                tmpFunc = str2func(['@' tmpIN sig.body tmpIN]);
                                [varargout{1:nargout}] = tmpFunc(in{:});
                            else
                                [varargout{1:nargout}] = this.func(in{:});
                            end
                        end
                    else
                        in = subs.subs;
                        computeUUID = uuidgen();
                        oPath = [formalFunc.storeLocation this.uuid filesep];
                        oFile = [computeUUID '_in.mat'];
                        fullFile = [oPath oFile];
                        mmkdir(oPath);
                        save(fullFile,'in');
                        inFile = argfile(fullFile,numel(in));
                    end
                case '{}'
                    
                    
                    %{
                    for e = 1:this.signature.size(2)
                        if e <= numel(subs(1).subs)
                            if ~isa(subs(1).subs{e},'fa')
                                X(e) = z(subs(1).subs{e});
                            else
                                X(e) = subs(1).subs;
                            end
                        else
                            X(e) = z0;
                        end
                    end
                    
                    
                    varargout{1} = m95(this.func,X);
                    %}
                    
                    
                    
                    %{
                    % forgot what this is for
                    % 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % get signature
                    sig = this.signature;
                    in = subs.subs;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % scan for classes
                    for e = 1:numel(in)
                        argClass{e} = class(in{e});
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % attach inputs which are portableQ's to the inputs of
                    % the function
                    stringConnect = find(strcmp(argClass,'portableQ'));
                    for e = 1:numel(stringConnect)
                        in{e}.addTarget(this);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % generate the output tokens
                    for e = 1:nargout
                        varargout{e} = portableQ();
                        varargout{e}.attachSource(this);
                    end
                    %}
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %queueVec = find(strcmp(argClass,'parallel.pool.DataQueue'));
                    
                    
                    %newIN = formalFunc.selectInputs(sig.pinIN,queueVec);
                    %gIN = formalFunc.buildGenericSignature('in','{}',numel(in));
                    %[newSignature] = formalFunc.replaceInputs(gIN,newIN(2:end-1),1);
                    
                    %{
                    queueVec = find(strcmp(argClass,'parallel.pool.DataQueue'));
                    
                    newIN = formalFunc.selectInputs(sig.pinIN,queueVec);
                    gIN = formalFunc.buildGenericSignature('in','{}',numel(in));
                    [newSignature] = formalFunc.replaceInputs(gIN,newIN(2:end-1),1);
                    tmpFunc = eval(['@' newIN sig.body newSignature]);
                    
                    % change the func call
                    this.func = tmpFunc;
                    this.updateSignature();
                    
                    % attach listener
                    afterEach(in{queueVec}, this.func);
                    %}
                    
            end
        end
        
        function [] = in_xform(this,ndim)
            
        end
        
        function [] = incrementVersion(this)
            this.version = this.version + 1;
        end
        
        function [] = initProfileTable(this)
            ts = datetime(datestr(now,'yyyy-mm-dd-HH-MM-ss'),'InputFormat','yyyy-MM-dd-HH-mm-ss');
           
            this.outputProfileTable = timetable(ts,NaN,NaN,NaN,NaN);
            this.outputProfileTable.Properties.VariableNames = ...
                {'Number_Returned_files','Expected_Time','Expected_Max_Memory','Max_Memory'};
            this.outputProfileTable(1,:) = [];
        end
        
        function [name] = getName(this)
            name = this.name;
        end
        
        function [body] = getFunctionBody(this)
            body = func2str(this.func);
            fidx1 = strfind(body,'(');
            fidx2 = strfind(body,')');
            while ~isempty(fidx1)
                body((fidx1(1)):(fidx2(1))) = [];
                fidx1 = strfind(body,'(');
                fidx2 = strfind(body,')');
            end
            body(1) = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % "OLD" way to eval, when subs ref was used for indexing stack
        function [result] = operate(this,varargin)
            result = this.func(varargin{:});
        end
        
        function [] = profile(this,collection,N)
            if nargin == 2;N = 10;end
            FileList = collection.getRandomSet(N);
            FileList = FileList.FileList;
            % loop over each file
            for e = 1:numel(FileList)
                fprintf(['start: dataFile:' num2str(e) ':' num2str(numel(FileList)) '\n'])
                % make temp location
                tmpLocation{e} = [tempname filesep];
                mmkdir(tmpLocation{e});
                
                % start clock and run
                strTime = clock;
                this.func(FileList{e},tmpLocation{e});
                stpTime = clock;
                
                % calculate delta time
                runTime(e) = etime(stpTime,strTime);
                
                % count return files
                cdir = dir(tmpLocation{e});cdir(1:2) = [];
                N(e) = numel(cdir);
                for f = 1:numel(cdir)
                    [~,~,extList{e}{f}] = fileparts(cdir(f).name);
                end
                fprintf(['end: dataFile:' num2str(e) ':' num2str(numel(FileList)) ':' num2str(runTime(e)) '\n'])
            end
            
            ts = datetime(datestr(now,'yyyy-mm-dd-HH-MM-ss'),...
                'InputFormat','yyyy-MM-dd-HH-mm-ss');
            this.outputProfileTable.Number_Returned_files(ts) = mode(N);
            this.outputProfileTable.Expected_Time(ts) = mean(runTime);
           
            
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % memory profile start
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            commenter = autoCommenter();
            [sig] = signature(this);
            [profileFile,profileFunction] = autoCommenter.copyFunctionFile(sig.body,'_PROFILE');
           
            [profileFuncHandle,profileFunctionName,profileFunctionFile] = ...
                commenter.moniterVars(profileFunction);
            
            
            % loop over each file
            for e = 1:numel(FileList)
                
                
                % make temp location
                tmpLocation{e} = [tempname filesep];
                mmkdir(tmpLocation{e});
                
                
                global totalMem;
                global totalCount;
                global lineN;
    
                totalMem = [];
                totalCount = [];
                lineN = [];
                
                profileFuncHandle(FileList{e},tmpLocation{e});
                
                
                maxMem(e) = max(totalMem);
            
                
                clear totalMem;
                clear totalCount;
                clear lineN;
            end  
            
            this.outputProfileTable.Expected_Max_Memory(ts) = mean(maxMem);
            this.outputProfileTable.Max_Memory(ts) = max(maxMem);
            %}
            
            
        end
        
        function [reportString,reportStruct] = lastProfile(this)
            preFIX = {'B','KB','MB','GB','TB','PB','EB'};
            
            if ~isempty(this.outputProfileTable)
            reportString = ['[' datestr(this.outputProfileTable.ts(end)) ' <---> ' ...
                num2str(this.outputProfileTable.Expected_Time(end))  ' <---> ' ...
                num2str(this.outputProfileTable.Number_Returned_files(end)) ']'];
            reportStruct.date = ['Last Profile Date: ' ...
                                datestr(this.outputProfileTable.ts(end))];
            reportStruct.expectedRunTime  = ['Expected Run Time: ' ...
                                num2str(this.outputProfileTable.Expected_Time(end))];
            reportStruct.expectedFilesReturned = ['Expected Number Output Files: ' ...
                                num2str(this.outputProfileTable.Number_Returned_files(end))];
                            
            memBytes = this.outputProfileTable.Expected_Max_Memory(end);
            if memBytes ~= 0
                ind = floor(log2(memBytes)/log2(1024));
                F = "%.2f";
                memBytes = compose(F,memBytes/1024^ind);
                reportStruct.expectedMaxMemory = ['Expected Max Memory: ' ...
                    char(memBytes)  preFIX{ind}];

                memBytes = this.outputProfileTable.Max_Memory(end);
                ind = floor(log2(memBytes)/log2(1024));
                F = "%.2f";
                memBytes = compose(F,memBytes/1024^ind);
                reportStruct.maxMemory = ['Max Memory: ' ...
                    char(memBytes) preFIX{ind}];
            else
                 reportStruct.expectedMaxMemory = ['Expected Max Memory: unknown'];
                 reportStruct.expectedMaxMemory = ['Max Memory: unknown'];
            end
            
            else
                reportString = ['[' 'NaN' ' <---> ' ...
                'NaN'  ' <---> ' ...
                'NaN' ']'];
                reportStruct.date = 'Last Profile Date: Never';
                reportStruct.expectedRunTime = 'Expected Run Time: Unknown';
                reportStruct.expectedFilesReturned = 'Expected Number Output Files: Unknown';
                reportStruct.expectedMaxMemory = 'Expected Max Memory: Unknown';
                reportStruct.maxMemory = 'Max Memory: Unknown';
            end
        end
        
        function [signature] = getSignature(this)
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % find the func
            signature = this.updateSignature();
        end
        
        function [signature] = updateSignature(this)
            try
                signatureString = func2str(this.func);
                sidx = strfind(signatureString,'@');
                eidx = strfind(signatureString,')');
                inputString = signatureString((sidx+1):(eidx(1)));
                inputString(1) = ',';
                inputString(end) = ',';
                n = numel(strfind(inputString,','))-1;
                inputString(1) = '(';
                inputString(end) = ')';
                this.signature.pinIN = inputString;
                this.signature.string = [this.signature.pinOUT '<-' this.signature.pinIN ']'];
                this.signature.size(2) = n;
            catch
                resetSignature(this)
            end
            signature = this.signature;
        end
        
        function [] = resetSignature(this)
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % find the func
            body = this.getFunctionBody();
            this.signature = formalFunc.signatureFromName(body);
        end
        
        
    end
    
    
    methods (Static)
    
        
        function [depList] = getDependencies(functionName)
            % find default function
            funcPath = which(functionName);
            % find list of required functions
            [fList,pList] = matlab.codetools.requiredFilesAndProducts(funcPath);
            % remove recursive calls and add to list
            if ~any(strcmp(fList,funcPath))
                fList{end+1} = funcPath;
            end
            % convert to file
            for e = 1:numel(fList)
                depList(e) = file(fList{e});
            end
        end
        
        
        function [signature] = signatureFromName(functionName)
            % find the .m file
            functionFile = which(functionName);
            % open the file
            fileID = fopen(functionFile);
            flag = true;
            % find the function line
            while flag
                line = fgets(fileID);
                line = strtrim(line);
                if contains(line,'function')
                    flag = false;
                end
            end
            fclose(fileID);
            % remove function word
            functionString = strrep(line,'function','');
            % trim the string
            functionString = strtrim(functionString);
            
            
            [signature] = formalFunc.signatureFromString(functionString);
            
        end
        
        function [signature] = signatureFromHandle(funcHandle)
            functionString = func2str(funcHandle);
            fidx0 = strfind(functionString,')');
            fidx1 = strfind(functionString,'(');
            functionBody = functionString((fidx0(1)+1):(fidx1(2)-1));
            signature = formalFunc.signatureFromName(functionBody);
        end
        
        function [signature] = signatureFromString(functionString)
            
            %
            eidx = strfind(functionString,'=');
            oidx = strfind(functionString,'(');
            functionName = functionString((eidx(1)+1):(oidx(1)-1));
            
            % remove body
            line = strrep(functionString,functionName,'');
            % output
            o1 = strfind(line,'[');
            o2 = strfind(line,']');
            out = line(o1:o2);
            % input
            i1 = strfind(line,'(');
            i2 = strfind(line,')');
            in = line(i1:i2);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            out = strrep(out,' ',',');
            in = strrep(in,' ',',');
            %%%%%%%%%%%%%%%%%%%%%%%%%
            outN = strfind(out,',');
            inN = strfind(in,',');
            signature.pinIN = in;
            signature.pinOUT = out;
            signature.string = [out '<-' in];
            signature.size = [numel(outN)+1 numel(inN)+1];
            signature.body = functionName;
        end
        
        function [func] = funcFromName(functionName)
            signature = formalFunc.signatureFromName(functionName);
            func = str2func(['@' signature.pinIN signature.body signature.pinIN]);
        end
        
        
        function [stringSignature] = removeEndInputs(stringSignature,n)
            stringSignature(1) = ',';
            stringSignature(end) = ',';
            fidx = strfind(stringSignature,',');
            stringSignature((fidx(n+1)+1):end) = [];
            stringSignature(1) = '(';
            stringSignature(end) = ')';
        end
        
        function [newSignature] = selectInputs(stringSignature,nvec)
            stringSignature(1) = ',';
            stringSignature(end) = ',';
            fidx = strfind(stringSignature,',');
            N = subs(fidx)-1;
            newSignature = [];
            for e = 1:numel(nvec)
                newSignature = [newSignature stringSignature(fidx(nvec(e)):fidx(nvec(e)+1))];
            end
            newSignature(1) = '(';
            newSignature(end) = ')';
        end
        
        function [newSignature] = replaceInputs(stringSignature,newInput,idx)
            stringSignature(1) = ',';
            stringSignature(end) = ',';
            fidx = strfind(stringSignature,',');
            strString = stringSignature(1:fidx(idx));
            stpString = stringSignature(fidx(idx+1):end);
            newSignature = [strString newInput stpString];
            newSignature(1) = '(';
            newSignature(end) = ')';
        end
        
        function [stringSignature] = buildGenericSignature(body,type,n)
            stringSignature = ',';
            for e = 1:n
               tmp = body;
               switch type
                   case '{}'
                       tmp = [tmp '{' num2str(e) '}'];
                   case '()'
                        tmp = [tmp '(' num2str(e) ')'];
                   case '_'
                        tmp = [tmp '_' num2str(e)];
                   case ''
                        tmp = [tmp num2str(e)];
               end
               stringSignature = [stringSignature tmp ','];
            end
            stringSignature(1) = '(';
            stringSignature(end) = ')';
        end
        
        
    end
    
end