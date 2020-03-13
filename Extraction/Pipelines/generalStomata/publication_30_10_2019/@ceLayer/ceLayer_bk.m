classdef ceLayer < generalizedComputeLayer

    
    properties
        %%%
        mode;
        %collectPtrPair;


        %%%
        ceType;
        pairLink;

        %%%
        upList;
        downList;
        leftList;
        rightList;
    

        %%%
        % collector store
        store;

    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % name := name of layer
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
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set pair link
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = setPairLink(obj,linkTo)
            obj.pairLink = linkTo;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % not yet needed
        function [] = compute(obj,varargin)
            fprintf('Not yet needed');
        end
        



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow logical for the collector/emitter 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flowData := data history, path commands (stop), etc
        function [result,flowData] = flow(obj,dFluid)

           


            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MIGHT RETURN HERE IF NOT WORK
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ret = obj.ceType;
            if ret == 0 & ~flowData.isTraceRoute
                flowData.traceRoute(obj,obj.downList{1}.bottomCE,'f')
                
                %
            end
            %}
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % report flow
            fprintf(['Calling flow on celayer:' obj.name '-' obj.uid '\n']);
            % call log with sqlite
            obj.logLayerLocation();
            % log 
            dFluid.logLocation(obj);
            % persist data flow
            obj.persistData(dFluid);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            switch dFluid.flowDirection
                case 'f'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % if flowing in forward direction
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %fprintf(['Flowing forward.']);


                    %%%%%%%%%%%%%%%%%%%%%%%%
                    % check for stopping conditions
                    % 1) at target
                    goFlag = ~dFluid.atCurrentTarget(obj);
                    % at edge of network
                    goFlag = goFlag & ~isempty(obj.downList);
                    %%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % check if this is a bottom node
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if strcmp(obj.ceType,'bottom')


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if in collect mode
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if strcmp(obj.mode,'collect')


                            
                            %%%%%%%%%%%%%%%%%%%%%%%%
                            % set the height/width ptr
                            %heightPtr = obj.collectPtrPair(1);
                            %widthPtr = obj.collectPtrPair(2);
                            % store the data
                            % OLD LINE BEFORE digital fluid was created and
                            % a handle object
                            %obj.store{heightPtr}{widthPtr} = dFluid;
                            % new line here
                            %obj.store{heightPtr}{widthPtr} = dFluid;
                            %here = 1 ;
                            %%%%%%%%%%%%%%%%%%%%%%%%
                            
                            result = {};
                            return;

                
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if in flow mode
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        elseif strcmp(obj.mode,'flow')




                            % if flow should continue - flow on
                            if goFlag



                                obj.emitterMode(dFluid);
                                %{
                                if size(obj.store,1) == 1
                                    [dFluid] = obj.downList{1}.flow(dFluid);
                                else
                                    [dFluid] = obj.downList{1}.flow(dFluid);
                                end
                                %}

                                %{
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % try auto fan  - width only
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                Width = numel(obj.downList);
                                for w = 1:Width
                                    [result,flowData] = obj.downList{w}.flow(flowData,varargin{:});
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %}


                             else

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % stop flow of computation and bounce back data
                                % bounce off the edge back to caller
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %result = varargin;
                                if size(obj.store,1) == 1
                                    if size(obj.store{1},2) == 1
                                        result = {obj.store{1}{1}{:}};
                                    else
                                        result = {obj.store{1}{:}};
                                    end
                                else
                                    result = obj.store;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            end

                        end


            
                    elseif strcmp(obj.ceType,'top')


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if go flag
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if goFlag
                            %{
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % try auto fan 
                            % each ceLayer has default double for loop
                            % the tall direction is the number of input args
                            % the wide direction is the number of layers
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % set collector pair to 'collect' mode
                            obj.pairLink.mode = 'collect';
                            inputStackHeight = size(dFluid,1);
                            Width = numel(obj.downList);
                            for h = 1:inputStackHeight
                                dFluid.splitStream(h,Width);
                                for w = 1:Width
                                    %%%%%%%%%%%%%%%%
                                    % ensure that after each callout
                                    % this layer is set to current
                                    obj.setAsCurrentLayer();


                                    %%%%%%%%%%%%%%%%
                                    dFluid.ptr = [h w];
                                    %tmpData = dFluid{h}{w};
                                    %tmpFluid = digitalFluid('tmp',tmpData);
                                    % set the pointer
                                    %obj.pairLink.collectPtrPair = [h w];
                                    % if inputStackHeight is one then only
                                    % single input call
                                    [dFluid] = obj.downList{w}.flow(dFluid);
                                end
                                % merge data from parallel streams
                                here = 1;
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %}

                            obj.emitterMode(dFluid);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % set the collector-pair into 'flow' mode
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %obj.pairLink.mode = 'flow';

                            % input to collector does not matter as it has
                            % collected the computation and will pass it on
                            [dFluid] = obj.pairLink.flow(dFluid);
                            
                        else

                            result = varargin;


                        end
                               

                    end


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'r'
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



        function [] = emitterMode(obj,dFluid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % try auto fan 
            % each ceLayer has default double for loop
            % the tall direction is the number of input args
            % the wide direction is the number of layers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set collector pair to 'collect' mode
            obj.pairLink.mode = 'collect';
            % get the height of the data
            inputStackHeight = size(dFluid,1);
            % get the width of the data
            Width = numel(obj.downList);
            % for each input
            for h = 1:inputStackHeight
                % create backup copy of fluid parcel
                bkFluidData = dFluid{h};
                %dFluid.splitStream(h,Width);
                for w = 1:Width
                    %%%%%%%%%%%%%%%%
                    % ensure that after each callout
                    % this layer is set to current
                    obj.setAsCurrentLayer();


                    %%%%%%%%%%%%%%%%
                    dFluid.ptr = [h w];
                    %tmpData = dFluid{h}{w};
                    %tmpFluid = digitalFluid('tmp',tmpData);
                    % set the pointer
                    %obj.pairLink.collectPtrPair = [h w];
                    % if inputStackHeight is one then only
                    % single input call
                    [dFluid] = obj.downList{w}.flow(dFluid);
                end
                % merge data from parallel streams
                here = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % restore ce-pair to flow
            obj.pairLink.mode = 'flow';
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % persist data function for calls to the CE layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = persistData(obj,flowData,varargin)
            %fprintf(['Calling data persist.\n']);
        
            if (flowData.persistData.toPersist == true)
                fprintf(['Data persist has been requested @' flowData.persistData.location '.\n'])
                switch flowData.persistData.location
                    case 'irods'
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make remote location location
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fprintf(['Making remote location.\n'])
                        rPath = ['/iplant/home/' flowData.persistData.userName ...
                                '/phytoLayers/className_' class(obj) ...
                                '/layerId_' obj.uid '/jobId_' flowData.jobData.uid '/'];
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



                end
            end
        end

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ceType
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ret] = ceType(obj)
                % 0 = top ---- 1 = bottom

            ret = 'None';

            if ~isempty(obj.downList) && ~isempty(obj.upList)
                ret = ~(~isa(obj.downList{1},'ceLayer') & isa(obj.upList{1},'ceLayer'));
            end

            if isempty(obj.upList)
                ret = false;
            end
        
            if isempty(obj.downList)
                ret = true;
            end

            if (ret == false)
                ret = 'top';
            else
                ret = 'bottom';
            end

        end
        %}
    end
    
end