classdef ceLayer <  uidTrackable & flowable & linkable
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
    %%%
    % ceLayer is now statless and memoryless. The states prior to Nov, 20
    % 2019 were flow and collect and the memory was the store for
    % collection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties
        % custom banner message
        bannerMessage;
        
        
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
        % collector store properties
        %mode;
        %store;          % store for data when in collector mode
        %storePointer;
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
            obj = obj@uidTrackable(name);
            % call the linkable class super constructor
            obj = obj@linkable();
            % set the type
            obj.ceType = type;


       

            %%%%%%%%%%%%%%%%%%%%
            % i think this can be removed per statless network 
            % set default mode to 
            % Removed Nov 26, 2019
            %obj.mode = 'flow';
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

                loggerText = ['flow-->' dFluid.flowDirection];
                
                
                dFluid = obj.logger(dFluid,loggerText);

                %%%%%%%%%%
                % check if persist is needed
                %%%%%%%%%%
                if dFluid.isGreen
                 %   dFluid.trigger(obj,1);
                end
                dFluid.trigger(obj,1);

            %{
            %%%%%%%%%%
            % logic for unwrap vs persist
            if dFluid.isBlue()
                dFluid = dFluid.trigger(obj,3);
            elseif dFluid.isRed()
                dFluid.trigger(obj,1);
            elseif dFluid.isGreen()
                dFluid.trigger(obj,1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % check for stopping conditions
                % 1) at target
                goFlag = ~any(dFluid.toStop(obj));
                % old syntax - REMOVE LATER
                %goFlag = ~dFluid.atCurrentStopTarget(obj);
                % at edge of network
                goFlag = goFlag & ~isempty(obj.downList);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the fluid color
                fluidColor = dFluid.color;
                % get the fluid direction
                fluidDirection = dFluid.flowDirection;
                % get the ceType
                ceType = obj.ceType;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%




                   

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % case fluid direction
                switch fluidDirection
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if forward
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'f'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % case ceType
                    switch ceType
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if bottomCE - forward
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'bottom'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if continue / stop
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if goFlag
                            if dFluid.isRed
                             %    dFluid = dFluid.trigger(obj,1);
                            end
                            % go
                            dFluid = obj.downList{1}.flow(dFluid);
                        else
                            % if there if more network to go
                            % update the fluid's stop queue
                            %if ~isempty(obj.downList)
                                % stop
                               dFluid.popStopTarget();
                            %end
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if topCE
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'top'
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if contine
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if goFlag
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % case fluid color
                            switch fluidColor
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % if green - forward
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            case "green"
                                dFluid = obj.emitterMode(dFluid);
                            case "red"
                                dFluid = dFluid.trigger(obj,1);
                                % check if fluid color is changed to green
                                if dFluid.isGreen()
                                    dFluid = obj.emitterMode(dFluid);
                                else
                                    % skip emitter - push on to next node
                                    dFluid = obj.downList{1}.flow(dFluid);
                                end
                            case "blue"
                                dFluid = dFluid.trigger(obj,3);
                                if dFluid.isWide
                                    dFluid = obj.emitterMode(dFluid);
                                else
                                    % skip emitter - push on to next node
                                    dFluid = obj.downList{1}.flow(dFluid);
                                end
                            end
                        end
                    end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if reverse
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'r'
                    [dFluid] = obj.reverseFlow(dFluid);
                end


            catch ME
                dFluid.errorLog = getReport(ME);
            end
            
        end

        function [] = printBanner(obj)
            
            
            
            if ~isempty(obj.bannerMessage)
                bannerLength = 80;
                
                
                spaceN = 4;
                nLevel = timingBlock('levelQuery');
                spaceBlock = repmat('   ',[1 nLevel]);
                bannerTop = repmat('*',[1 bannerLength - numel(spaceBlock)]);
                bannerTop(1:2) = [];
                bannerTop = ['|' bannerTop '|'];
                    
                bannerBuffer = repmat(' ',[1 bannerLength - numel(spaceBlock)]);
                bannerBuffer(1:4) = [];
                bannerBuffer = ['|*' bannerBuffer '*|'];

                fprintf([bannerTop '\n']);
                fprintf([bannerBuffer '\n']);
                
                for block = 1:numel(obj.bannerMessage)

                   
                    
                    bannerMessage = ['|*    ' spaceBlock obj.bannerMessage{block} ' '];




                    

                    while numel(bannerMessage) > 0

                        fidx = strfind(bannerMessage,' ');
                        lineLen = (bannerLength-2*spaceN+2);
                        lidx = find((fidx - lineLen) <= 0);
                        sidx = fidx(lidx(end));




                        tmpMsg = bannerMessage(1:sidx);
                        bannerMessage(1:sidx) = [];
                        if ~isempty(bannerMessage)
                            bannerMessage = ['|*    ' spaceBlock bannerMessage];
                        end
                        tmpSpace = repmat(' ',[1 bannerLength - numel(tmpMsg)]);
                        tmpSpace(1:2) = [];

                        tmpMsg = [tmpMsg tmpSpace '*|'];


                        fprintf([tmpMsg '\n']);
                    end
                    
                    
                end
                fprintf([bannerBuffer '\n']);
                fprintf([bannerTop '\n']);
                
                %bannerTop = repmat('*',[1 bannerLength - numel(spaceBlock)]);
                %fprintf([spaceBlock bannerTop]);
                
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % logger, persister, unwrapper
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = logger(obj,dFluid,functionCall)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call superclass logger
            logger@uidTrackable(obj,functionCall);
            
            
            % report flow
            ceType = obj.ceType;
            msg = ['[' obj.name '@' obj.uid ']'];
            switch ceType
                case 'bottom'
                    timingBlock('stop',msg);
                case 'top'
                    obj.printBanner();
                    timingBlock('l-start',msg);
            end
            
            
            % log with fluid
            dFluid.logger(obj);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % color change
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            % prepare the width of the input
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if dFluid.isWide()
                here = 1;
            end


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
                        if obj.isWide()
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
                                if obj.isWide()
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
                                if obj.isWide()
                                    obj.pairLink.mode = 'flow';
                                end

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                % call the bottom ce pair to continue
                                if obj.isWide()
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % not parallel
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for each input
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for h = 1:inputStackHeight
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if topCE is not simple then bottom should
                    % also be complex
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if obj.isWide()
                        for w = 1:Width

                            if ~dFluid.isWide()
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % make working fluid from backup
                                % new parcel for each branch
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                wFluid(w) = copy(dFluid);
                                wFluid(w).pushStopTarget(obj.pairLink);
                            else
                                wFluid(w) = dFluid(w);
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % call thread until bottom ce-pair is reached
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            wFluid(w) = obj.downList{w}.flow(wFluid(w));
                        end

                        dFluid = obj.fuse(wFluid);
                        [dFluid] = obj.pairLink.flow(dFluid);


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % topCE is not wide
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else

                        dFluid = obj.downList{1}.flow(dFluid);


                        here = 1;



                    end

                end
            end


           
        end


        function [dFluid] = prepareFluidWidth(dFluid)

        end


    end
    


    methods (Access=private)


        function [b] = isWide(obj)
            if strcmp(obj.ceType,'bottom')
                b = numel(obj.upList) ~= 1;
            elseif strcmp(obj.ceType,'top')
                b = numel(obj.downList) ~= 1;
            end
        end

        %{
        function [] = collect(obj,newData)
            %fprintf(['Calling collect mode.Layer name:' obj.name '.\n']);
            %newData{1}
            obj.store{obj.storePointer} = obj.glue(obj.store{obj.storePointer},newData);
        end
        %}


        function [dFluid] = fuse(obj,wFluid)
            if wFluid.isGreen()

                for e = 1:(numel(wFluid)-1)
                    wFluid(e+1).data = obj.glue(wFluid(e).data,wFluid(e+1).data);
                end
                dFluid = wFluid(end);

            else

                for e = 1:numel(wFluid)
                    if wFluid(e).isGreen()
                        wFluid(e).toCompute = false;
                    end
                end
                % make complex
                dFluid = digitalFluid('blueWrap',wFluid,"blue");
                dFluid.assignUnwrapNode(obj);
                
            end
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