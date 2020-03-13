classdef persistQueue < eventQueue

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % persistQueue constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = persistQueue()
            obj@eventQueue('persistQueue',false);
        end
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % event function - will save the digital fluid 
        % local or remote 
        % flips the mode to read and the color to red
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = eventFunction(obj,curNode,eventIdx,dFluid)
            fprintf(['Data persist has been requested @' obj.eventData(eventIdx).location '.\n'])
            switch obj.eventData(eventIdx).location
                case 'irods'
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % make remote location location
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    fprintf(['Making remote location.\n'])
                    rPath = ['/iplant/home/' obj.eventData(eventIdx).userName ...
                            '/phytoLayers/className_' class(obj) ...
                            '/layerId_' obj.uid '/jobId_' obj.jobData.uid '/'];
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

                    if strcmp(obj.eventData(eventIdx).ioType,'write')
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make local location
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fprintf(['Making local location.\n'])
                        localBase = ['~' filesep 'layerCache' filesep];
                        lPath = [localBase ...
                                'phytoLayers/className_' class(obj) ...
                                '/layerId_' curNode.uid '/jobId_' obj.eventData(eventIdx).jobID '/'];
                        CMD = ['mkdir -p ' lPath];
                        r = system(CMD);


                        [~,tmpName] = fileparts(tempname);
                        localFileName = [lPath tmpName '.mat'];
                        save(localFileName,'dFluid');


                        dFluid.setCurrent(localFileName);
                        dFluid.color = "red";

                        obj.eventData(eventIdx).ioType = 'read';


                    elseif strcmp(obj.eventData(eventIdx).ioType,'read')
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                        fileName = dFluid.getCurrent();
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%% use color conversion tool
                        
                        %%% try manual setting here of loaded object
                        loadedObj = load(fileName,'dFluid');
                        loadedObj = loadedObj.dFluid;


                
                        % set the data
                        dFluid.setCurrent(loadedObj.getCurrent());
                        % set the color to green
                        dFluid.color = "green";
                        % pop the event from the queue
                        obj.pop(eventIdx);
                        here = 1;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define trigger - call super class trigger
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dFluid] = trigger(obj,curNode,dFluid)
            dFluid = trigger@eventQueue(obj,curNode,dFluid);
        end

    end


    methods (Static)


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct the persist data struct
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [eventData] = generatePersistData(type,user,jobID)
            if isnumeric(jobID);jobID = num2str(jobID);end
            eventData.location = type;
            eventData.ioType = 'write';
            eventData.user = user;
            eventData.jobID = jobID;
        end



    end


end