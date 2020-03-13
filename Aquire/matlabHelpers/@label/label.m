classdef label < tile
    % a lable isa tile and associated human readable text
    properties
        hText; % humand readable text
    end
    
    methods
        % constructor for label
        function [obj] = label(hText)
            % make the tile
            obj = obj@tile();
            % attach the human text
            obj.hText = hText;
        end
    end
    
end
