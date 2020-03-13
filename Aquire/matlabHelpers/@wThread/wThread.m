classdef wThread < doid
    % worker thread - isa doid
    properties
        
    end
    
    
    methods
        
        
        function [obj] = wThread(sharedToMain)
            % shared channel to main
            obj.sharedToMain = sharedToMain;
            % shared in channel
            obj.sharedIn = parallel.pool.DataQueue;
            % attach this function as a listener from main
            afterEach(obj.sharedIn, @(X)obj.processMessage(X));
            % transmitt ID
            obj.sendID();
        end
        
        function [] = processMessage(obj,msg)
            % ensure the msg was meant for this
            %if msg.to == oid.proj(obj)
                %retMsg = etime(clock,msg);
                obj.halt();
            %end
        end
        
        
        function [] = sendSharedToMain(obj,msg)
            % send the channels back to main
            send(obj.sharedToMain,msg);
        end
        
        function [msg] = buildDataMessage(obj,to,data)
            % from this
            from = oid.proj(obj);
            % container
            body = containerMsg(data);
            % build msg
            msg = daqMessage(from,to,body,true);
        end
        
        function [] = sendDataMessage(obj,to,data)
            % build the message
            msg = obj.buildDataMessage(to,data);
            % send the messsage
            obj.sendSharedToMain(msg);
        end
        
        
        function [] = sendID(obj)
            % make addressBook entry and send it back to main
            addressEntry = addressBookEntry(obj.uuid,obj.sharedIn);
            % send the data back to the main
            obj.sendDataMessage('',addressEntry);
        end
        
    end
end