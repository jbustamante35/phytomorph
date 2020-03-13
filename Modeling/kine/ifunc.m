function [P,L] = ifunc(P,F,W,V,L)
    try
        % X = init point
        % F = function to eval
        % W = window to eval function - width and number of points
        % V = search function - somehow
        % L = zoom loop number
        if ~L.halt(L.totalLoop,L.resolutionValue)

        %for l = 1:L.totalLoop
            idx = [];
            while numel(idx) ~= 1
                % test call for simple zoom
                [Y,X] = zoomOnPoint(F,P,W);
                % search for the point 
                idx = V(Y);
                % increase the number of points if there are many maxs
                %if numel(idx) > 1;W.numP = 2*W.numP;end
                if numel(idx) > 1
                    return;
                end
            end
            % assign the new point
            P = X(idx,:);
            % update the search window
            W.width = W.width*L.loopCompression;
            % update the loop count
            L.totalLoop = L.totalLoop + 1;
            % update resolution
            if L.totalLoop > 1
             L.resolutionValue = L.updateCurrentResolution(P,L.history);
            end
            % update the history
            L.history = [L.history;P];

            % recurse
            [P,L] = ifunc(P,F,W,V,L);
        end
    catch ME
        ME
    end
end