classdef ceLayer < generalizedComputeLayer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ceLayer is a collector/emitter node.  Each device layer has two of
    % these, the input and output pair.  Each ceLayer can serve as a
    % collector and emitter.  If the node is in collector mode, then it
    % will use the store to collect the data passed in. This design element
    % might/will change soon.  The fluid will be 'upgraded' to have the
    % statefullness and logic rather than embedding statefullnes in the
    % network.  Making the network passive (non-statefull) implies that any
    % copy of the network can be used rather than saving a state copy of
    % the network.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties
        % the threshold for parallelization
        distributionThresh = 10;

       
        % funnel commands
        funnelCommands;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the type of ce-layer top or bottom
        ceType;
        % link to the top/bottom pair
        pairLink;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set of lists for pointers
        upList;
        downList;
        leftList;
        rightList;
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % collector store properties
        mode;
        store;          % store for data when in collector mode
        storePointer;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % name := name of layer
        % type := top/bottom type
        function [obj] = ceLayer(name,type)
            % call to super class
            obj = obj@generalizedComputeLayer(name);
            % list connections to layers
            obj.upList = {};
            obj.downList = {};
            obj.leftList = {};
            obj.rightList = {};
            obj.store = {};
            % set the type
            obj.ceType = type;
            % set default mode to 
            obj.mode = 'flow';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set pair link
        % each top and bottom are paired
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setPairLink(obj,linkTo)
            obj.pairLink = linkTo;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute
        % the ceLayer does not yet call or need compute
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % not yet needed
        function [] = compute(obj,varargin)
            fprintf('Not yet needed');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % commands to funnel/fuse data from parallel flows
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setFunnelCommands(obj,fCommands)
            obj.funnelCommands = fCommands;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow logical for the collector/emitter 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flowData := data history, path commands (stop), etc
        function [dFluid] = flow(obj,dFluid)
            try
                obj.logger(dFluid,['flow-->' obj.mode '-->' dFluid.flowDirection]);
                switch dFluid.flowDirection
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if flowing in forward direction
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'f'
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %fprintf(['Flowing forward.']);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % check for stopping conditions
                        % 1) at target
                        goFlag = ~dFluid.atCurrentTarget(obj);
                        % at edge of network
                        goFlag = goFlag & ~isempty(obj.downList);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % check if this is a bottom node
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if strcmp(obj.ceType,'bottom')
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % if in collect mode
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if strcmp(obj.mode,'collect')
                                obj.collect(dFluid{dFluid.ptr})
                                %obj.store{end+1} = dFluid{1};
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % if in flow mode
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            elseif strcmp(obj.mode,'flow')
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % if collector is complex
                                % complex means fusing more than one data
                                % stream
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if ~obj.isWide()
                                    dFluid = obj.hydrateStore(dFluid);
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % if flow should continue - flow on
                                % design elements demand that downlist from
                                % bottom ceLayer can only have one
                                % connection
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if goFlag
                                    dFluid = obj.downList{1}.flow(dFluid);
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if go flag
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        elseif strcmp(obj.ceType,'top')
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % if contine
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if goFlag
                                dFluid = obj.emitterMode(dFluid);
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if flowing in reverse direction
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'r'
                        [dFluid] = obj.reverseFlow(dFluid);
                end
            catch ME
                report = getReport(ME);
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = logger(obj,dFluid,functionCall)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call superclass logger
            logger@generalizedComputeLayer(obj,functionCall);
            % log with fluid
            dFluid.logger(obj);
            % persist data flow
            obj.persistData(dFluid);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % not needed and will be cleared 
        function [] = clearStore(obj)
            obj.store{obj.storePointer} = {};
        end

        function [wFluid] = hydrateStore(obj,wFluid)
            wFluid.data{obj.storePointer} = obj.store{obj.storePointer};
            obj.clearStore();
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sequence to stream data 
        % from a emitter,through computation layer, 
        % to collector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = emitterMode(obj,dFluid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % try auto fan 
            % each ceLayer has default double for loop
            % the tall direction is the number of input args
            % the wide direction is the number of layers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the height of the data
            inputStackHeight = size(dFluid,1);
            % get the width of the data
            Width = numel(obj.downList);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if simple ceLayer - then call flow-compute
            % else set to collect first
            if ~obj.isWide()
                % set collector pair to 'collect' mode
                obj.pairLink.mode = 'collect';
                % set the collector store to empty cell array
                obj.pairLink.store = cell(inputStackHeight,1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if inputStackHeight > (obj.distributionThresh - 1) && dFluid.canParallel
                %fprintf(['Calling parallel process.Layer Name:'  obj.name '\n']);
            
                %parpool(3)

                

                
               
                spmd
                    
                    codist = codistributor1d(1,[],[inputStackHeight,1]);
                 
                    labindex;
    
                    inputArray = codistributed(dFluid.data,codist);
                    localData = getLocalPart(inputArray);
                    


                    localData;

                    pFluid = digitalFluid('parallelFluid',localData);
                    pFluid.stopData.targetQueue = dFluid.stopData.targetQueue;
                    pFluid.toCompute = dFluid.toCompute;
                
                    LocalinputStackHeight = size(pFluid);

                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % for each input
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for h = 1:LocalinputStackHeight
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if topCE is not simple then bottom should
                        % also be complex
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if ~obj.isWide()
                            %fprintf(['Calling Wide branch.\n']);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % set collector for store pointer
                            obj.pairLink.storePointer = h;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % create backup copy of fluid parcel
                            bkFluidData = pFluid{h};



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLOW MODE NOW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if ~obj.isWide()
                                    % set collector pair to 'collect' mode
                                    obj.pairLink.mode = 'collect';
                                end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                            for w = 1:Width




                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % ensure that after each callout
                                % this layer is set to current
                                % for tracing network
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %obj.setAsCurrentLayer();


                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % make working fluid from backup
                                % new parcel for each branch
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                wFluid = digitalFluid('workingFluid',{bkFluidData});
                                wFluid.stopData.targetQueue = pFluid.stopData.targetQueue;
                                wFluid.toCompute = pFluid.toCompute;

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % call thread until bottom ce-pair is reached
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                wFluid = obj.downList{w}.flow(wFluid);




                            end





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLOW MODE NOW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % restore ce-pair to flow if complex
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if ~obj.isWide()
                                    obj.pairLink.mode = 'flow';
                                end

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                % call the bottom ce pair to continue
                                if ~obj.isWide()
                                    [pFluid] = obj.pairLink.flow(pFluid);
                                end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FLOW MODE NOW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % topCE is not wide
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % create backup copy of fluid parcel
                            %bkFluidData = dFluid{h};

                            
                            wFluid = digitalFluid('workingFluid',{pFluid{h}});
                            wFluid.stopData.targetQueue = pFluid.stopData.targetQueue;
                            wFluid.toCompute = pFluid.toCompute;
                            wFluid = obj.downList{1}.flow(wFluid);
                            pFluid{h} = wFluid{1};
                            
                            %{
                            
                            pFluid.ptr = h;
                            pFluid = obj.downList{1}.flow(pFluid);
                            %}
                        end

                        pFluid{1};
                        

                    end
                    

                    here = 1;


                    localData = pFluid.data;
                    %done = gather(localData,codist);
                    done = codistributed.build(localData,codist);

                end

                dFluid.data = gather(done,codist);


            else
            

                
                %{
                fprintf(['Calling serial process.Layer Name:'  obj.name '\n']);
                if iscell(dFluid{1}{1})
                    a = dFluid{1}{1};
                    a{1}
                else
                    a = dFluid{1}{1};
                    a(1)
                end
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for each input
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for h = 1:inputStackHeight
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if topCE is not simple then bottom should
                    % also be complex
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if ~obj.isWide()
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % create backup copy of fluid parcel
                        bkFluidData = dFluid{h};


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % set collector for store pointer
                        obj.pairLink.storePointer = h;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FLOW MODE NOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if ~obj.isWide()
                                % set collector pair to 'collect' mode
                                obj.pairLink.mode = 'collect';
                            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                        for w = 1:Width
                            %{
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % set collector for store pointer
                            obj.pairLink.storePointer = [h w];
                            %}

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % ensure that after each callout
                            % this layer is set to current
                            % for tracing network
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            obj.setAsCurrentLayer();


                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % make working fluid from backup
                            % new parcel for each branch
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            wFluid = digitalFluid('workingFluid',{bkFluidData});
                            wFluid.stopData.targetQueue = dFluid.stopData.targetQueue;
                            wFluid.toCompute = dFluid.toCompute;

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % call thread until bottom ce-pair is reached
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            wFluid = obj.downList{w}.flow(wFluid);




                        end





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FLOW MODE NOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % restore ce-pair to flow if complex
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if ~obj.isWide()
                                obj.pairLink.mode = 'flow';
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                            % call the bottom ce pair to continue
                            if ~obj.isWide()
                                [dFluid] = obj.pairLink.flow(dFluid);
                            end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FLOW MODE NOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % topCE is not wide
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % create backup copy of fluid parcel
                        %bkFluidData = dFluid{h};

                        %{
                        wFluid = digitalFluid('workingFluid',{dFluid{h}});
                        wFluid.stopData.targetQueue = dFluid.stopData.targetQueue;
                        wFluid.toCompute = dFluid.toCompute;
                        wFluid = obj.downList{1}.flow(wFluid);
                        dFluid{h} = wFluid{1};
                        %}


                        dFluid.ptr = h;
                        dFluid = obj.downList{1}.flow(dFluid);




                    end




                end







            end


           
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % restore ce-pair to flow if complex
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~obj.isWide()
                obj.pairLink.mode = 'flow';
            end
            %}

            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % if not simple then call the 
            % bottom ce pair to continue
            if ~obj.isWide()
                [dFluid] = obj.pairLink.flow(dFluid);
            end
            %}
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % persist data function for calls to the CE layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = persistData(obj,dFluid,varargin)
            %fprintf(['Calling data persist.\n']);
            toPersist = dFluid.persistData.toPersistFunction(obj.uid);
            if (toPersist == true)
                fprintf(['Data persist has been requested @' dFluid.persistData.location '.\n'])
                switch dFluid.persistData.location
                    case 'irods'
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make remote location location
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fprintf(['Making remote location.\n'])
                        rPath = ['/iplant/home/' dFluid.persistData.userName ...
                                '/phytoLayers/className_' class(obj) ...
                                '/layerId_' obj.uid '/jobId_' dFluid.jobData.uid '/'];
                        CMD = ['imkdir -p ' rPath];
                        r = system(CMD);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % saving temp file locaal
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fprintf(['Saving local data for temp pre transfer.\n']);
                        % default tmp location is:
                        tmpLocation = '~/phytoMorphTK/cache/';
                        CMD = ['mkdir -p ' tmpLocation];
                        system(CMD);
                        [~,tmpName] = fileparts(tempname);
                        localFileName = [tmpLocation filesep tmpName '.mat'];
                        remoteFileName = [rPath tmpName '.mat'];
                        save(localFileName,'varargin');
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % pushing data
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        CMD = ['iput -f ' localFileName ' ' remoteFileName];
                        system(CMD);
            
                    case 'local'
                        here = 1;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make remote location location
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fprintf(['Making remote location.\n'])
                        localBase = ['~' filesep 'layerCache' filesep];
                        lPath = [localBase ...
                                'phytoLayers/className_' class(obj) ...
                                '/layerId_' obj.uid '/jobId_' dFluid.jobData.uid '/'];
                        CMD = ['mkdir -p ' lPath];
                        r = system(CMD);


                        [~,tmpName] = fileparts(tempname);
                        localFileName = [lPath tmpName '.mat'];
                        save(localFileName,'dFluid');


                end
            end
        end

    end



    methods (Access=private)


        function [b] = isWide(obj)
            if strcmp(obj.ceType,'bottom')
                b = numel(obj.upList) == 1;
            elseif strcmp(obj.ceType,'top')
                b = numel(obj.downList) == 1;
            end
          
        end


        function [] = collect(obj,newData)
            %fprintf(['Calling collect mode.Layer name:' obj.name '.\n']);
            %newData{1}
            obj.store{obj.storePointer} = obj.glue(obj.store{obj.storePointer},newData);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % glue commands - default to glue as you go
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [out] = glue(obj,in1,in2)
            if isempty(in1)
                out = in2;
            else
                for e = 1:numel(in1)
                    %curType = obj.funnelCommands{e};
                    curType = ceLayer.buildGlueCommand(in1{e},in2{e});
                    if strcmp(curType{1},'cat')
                        out{e} = obj.catType(in1{e},in2{e},curType{2});
                    elseif strcmp(curType{1},'keep')
                        out{e} = obj.keepType(in1{e},in2{e},curType{2});
                    end
                end
            end
        end




        function [out] = catType(obj,in1,in2,dim)
            out= cat(dim,in1,in2);
        end

        function [out] = keepType(obj,in1,in2,select)
            if (select == 1)
                out = in1;
            elseif (select == 2)
                out = in2;
            end
        end



        function [dFluid] = reverseFlow(obj,dFluid)

    
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if flowing in reverse direction
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %fprintf(['Flowing reverse.']);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % check to see if flow should be reversed because the
                        % top/edge of the network has been hit
                        if isempty(obj.upList)
                            dFluid.flowDirection = 'f';
                            [dFluid] = obj.downList{1}.flow(dFluid);
                        else
                            [dFluid] = obj.upList{1}.flow(dFluid);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end

    end


    methods (Static)
        %{
        function [cmd] = defaultCommands(width)
            for e = 1:width
                cmd{e} = {'cat',1};
            end
        end
        %}

    
        function [cmd] = buildGlueCommand(in1,in2)
            if ~strcmp(class(in1),class(in2))
                cmd{1} = 'cat';
                cmd{2} = 2;
            else
                if ischar(in1)
                    cmd{1} = 'keep';
                    cmd{2} = 1;
                else
                    cmd{1} = 'cat';
                    cmd{2} = 2;
                end
            end
        end
    end


end