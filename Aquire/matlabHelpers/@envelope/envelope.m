classdef envelope < doid
    
    properties
        to;         % what/who the message is to 
        message;    % the message
    end
    
    methods
        % constructor for the envelope
        % note that the constructor needs to handle the 
        % zero input case for the datastore re-construction
        function [obj] = envelope(to,message)
            if nargin ~= 0
                % transfor the to property into pointer
                obj.to = dptr.transformInput(to);
                % transfor the to property into pointer
                obj.message = dptr.transformInput(message);
            end
        end
    end
    
end