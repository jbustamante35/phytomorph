classdef stateCube < messageRouter
    
    
    properties
        % type
        type;
        
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
        
       
        
        type,startFunction,channelGenFunc,
        
        parentAddressEntry
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = stateCube(parentAddressEntry,childN,pauseDuration)
            
            % generate the self channel
            inChannel = channelGenFunc{1}();
            % init the state cubde as a message router
            obj = obj@messageRouter(inChannel);
            % init the type
            obj.type = type(1);
            % add parent as first entry to address book
            obj.addEntry(parentAddressEntry);
            % create the pool if one does not exist
            if childN(1) ~= 0
                obj.pool = gcp('nocreate');
                if isempty(obj.pool)
                    obj.pool = parpool('local',childN(1));
                end
            end
            % grab the first start function
            if ~isempty(startFunction);sFunc = startFunction(1);end
            
            
            
            
            % pop the used child data
            type = type(2:end);
            startFunction = startFunction(2:end);
            channelGenFunc = channelGenFunc(2:end);
            selfEntry = obj.getSelfEntry();
            childN = childN(2:end);
            pauseDuration = pauseDuration(2:end);
            
            
            
            % eval par on pool
            future = parfevalOnAll(obj.pool,@(A,B,C,D)stateCube(A,B,C,D),1,...
                channelGenFunc,selfEntry,childN,pauseDuration);

            
            
            
             deltaDelay = 30; % set the delay to 3 sec
                publicFromWorkers = parallel.pool.DataQueue;
                
                
                
           
            % attach this function as a listener from main
            afterEach(publicFromWorkers, @(X)obj.processMessage(X));
            
            
            % incoming from all
            sharedIn;
            
            
            % set the halt flag to false
            obj.isAlive = true;
            
            
            % set the worker into an inf loop with a stop condition
            obj.go(pauseDuration);
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