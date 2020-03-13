classdef tile < doid
    % a tile is a doid with an envolope
    % the envelope contains the tile's message
    properties
        envelope; % envelope for the tile
    end
    
    methods
        % constructor
        function [obj] = tile()
            obj = obj@doid();
        end
        % attach message - should be in envelope
        function [] = attachEnvelope(obj,envelope)
            % message for tile - make pointer
            obj.envelope = dptr.transformInput(envelope);
        end
        % generate the image for the message via envelope id
        function [I] = qrImage(obj,sz)
            % the QR-image will contain the uuid of an envelope
            I = generateQRtile(obj.envelope.refs.uuid,sz);
            % it could contain the uuid of the pointer to the msg
            
        end
        % get the tile type
        function [type] = tileType(obj)
            type = obj.envelope.messageType();
        end
    end
end