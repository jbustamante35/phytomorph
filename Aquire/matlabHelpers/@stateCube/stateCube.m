classdef stateCube < messageRouter
    %%%%%%%%%%%%%%%%%%
    % question: May 11 2020
    % Can a paused thread still execute a function when the queue it
    % populated ?
    %%%%%%%%%%%%%%%%%%
    % As of May 11 2020 - the state cube will have a depth of one only.
    % This implies that a sub process can not have a subprocess.  The local
    % toplogy is a fan out only.  The main can connect to another process
    % and share data via channel.  Two processes currently need to
    % communicate through an intermediate process.  aka. communication is
    % hierichal and not a network.
    %%%%%%%%%%%%%%%%%%
    properties
        % type
        %type;
        
        % alive status
        isAlive;
        
        % router for communication
        router;
        
        % access to pool
        pool;
        
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % channelGenFunc := generate channel
        % parent address
        % pause duration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = stateCube(parentChannel,parentAddressEntry,childN,pauseDuration)
            
            % pass the constructor for the channel Generator function
            % to the stateCube
            
            
            
            % generate the self channel
            %inChannel = channelGenFunc();
            % init the state cube as a message router
            %obj = obj@messageRouter(inChannel);
            obj = obj@messageRouter('');
            
            
            
            % init the type
            %obj.type = type(1);
            
            % add parent as first entry to address book
            %obj.addEntry(parentAddressEntry);
            
            
            % create the pool if one does not exist
            if childN ~= 0
                
                obj.pool = gcp('nocreate');
                if isempty(obj.pool)
                    obj.pool = parpool('local',childN(1));
                end
                
                for e = 1:childN
                    backChannel{e} = parallel.pool.DataQueue;
                    
                    % attach this function as a listener from main
                    afterEach(backChannel{e}, @(X)obj.processMessage(X));
                    
                end
                
                NC = mat2cell(zeros(1,childN), [1],ones(1,childN));
                PD = mat2cell(pauseDuration*ones(size(backChannel)), [1],ones(1,childN));
            
                % eval par on pool
                for c = 1:childN
                    future{c} = parfevalOnAll(obj.pool,@(A,B,C,D)stateCube(A,B,C,D),1,...
                        backChannel{c},[],NC{c},PD{c});
                end
                
            end
            
            % grab the first start function
            %if ~isempty(startFunction);sFunc = startFunction(1);end
            
            
           
            %{
            % pop the used child data
            type = type(2:end);
            startFunction = startFunction(2:end);
            channelGenFunc = channelGenFunc(2:end);
            selfEntry = obj.getSelfEntry();
            childN = childN(2:end);
            pauseDuration = pauseDuration(2:end);
            %}
            

            
            
            
            %deltaDelay = 30; % set the delay to 3 sec
            
            %
            %publicFromWorkers = parallel.pool.DataQueue;
                
           
            % incoming from all
            %sharedIn;
            
            
            % set the halt flag to false
            obj.isAlive = true;
            
            
            % set the worker into an inf loop with a stop condition
            if ~isnan(pauseDuration);obj.go(pauseDuration);end
        end
        
        
        
        function [] = forceState(obj,state)
        
        end
        
        
        function [] = configure(obj,ID)
            obj.ID = ID;
        end
        
        
        function [] = halt(obj)
            obj.isAlive = false;
        end
        
        function [] = go(obj,pauseDuration)
            
           % send(obj.sharedToMain,obj.privateFromMain);
            
            while obj.isAlive
                pause(pauseDuration);
            end
            
            
        end
        
    end
end

%{

    obj = stateCube(0,'',4,30);
%}