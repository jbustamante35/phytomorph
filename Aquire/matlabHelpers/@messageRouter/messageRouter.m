classdef messageRouter < doid

    properties
        % address book
        addressBook;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor for a router
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = messageRouter(inChannel)
            % init object
            obj = obj@doid();
            % init address book - add self
            obj.addressBook = addressBook(obj.uuid,inChannel);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add entry to address book
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = addEntry(obj,entry)
            % add entry to address book
            obj.addressBook.addEntry(entry);
        end
        
        function [entry] = getSelfEntry(obj)
            entry = obj.addressBook.getSelfEntry();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = processMessage(obj,msg)
            here = 1;
            %
            if isa(class(message.body.data),'function_handle')
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % send message to parent
        % parent is #2 in list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = sendToParent(obj,msg)
            % send message back to main
            send(obj.addressBook(2).channel,msg);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build a data message
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [msg] = buildDataMessage(obj,to,data)
            % from this
            from = oid.proj(obj);
            % container message body
            body = containerMsg(data);
            % build msg
            msg = daqMessage(from,to,body,true);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % send data message
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = sendDataMessage(obj,to,data)
            % build the message
            msg = obj.buildDataMessage(to,data);
            % send the messsage
            obj.sendSharedToMain(msg);
        end
        
        
    end
    
    
    
end