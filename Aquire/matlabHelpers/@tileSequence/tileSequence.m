classdef tileSequence < tile
    properties
        title;
        metaN = 2; % start and stop
        sampleN = 1;  % default to one sample QR code
        
        usbPort;
        tileSeq;
    end
    
    methods
        
        
        function [obj] = tileSequence(title,usbPort)
            obj.title = title;
            if nargin > 1;obj.usbPort = dptr.transformInput(usbPort);end
            obj.tileSeq = oid.empty(0,3);
        end
        
        
        function [] = attachStartTile(obj,strMsg)
            obj.tileSeq(1) = strMsg;
        end
        
        
        function [] = attachEndTile(obj,endMsg)
            obj.tileSeq(2) = endMsg;
        end
        
        
        
        
        
        
        %{
        function [] = generate(obj,strBody)
            
            global store;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate start tile
            strTile = store.generate('tile');
            % generate message for start tile
            %strMessage = daqMessage(strTile,obj.usbPort,strMsg);
            strMessage = store.generate('daqMessage',strTile,obj.usbPort,strBody);
            % attach message for start tile
            store.invoke(strTile,'attachMessage',strMessage);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate stop tile
            %stpTile = tile();
            stpTile = store.generate('tile');
            % generate end msg
            endBody = endMsg();
            % generate message for stop tile
            %stpMessage = daqMessage(stpTile,obj.usbPort,endBody);
            stpMessage = store.generate('daqMessage',strTile,obj.usbPort,endBody);
            % attach message to stop tile
            store.invoke(strTile,'attachMessage',stpMessage);
            
            obj.tileSeq(1) = dptr(strTile);
            obj.tileSeq(2) = dptr(stpTile);
            
            %
            store.write(obj);
        end
        %}
    end
end