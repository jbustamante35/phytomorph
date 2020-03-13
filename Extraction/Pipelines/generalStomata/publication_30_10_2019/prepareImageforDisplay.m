function [I] = prepareImageforDisplay(I,prepareType,varargin)
    switch prepareType
        case 'edgeTrim'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for r = 1:4
            I(:,1:varargin(1)) = [];
            I = imrotate(I,90);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end