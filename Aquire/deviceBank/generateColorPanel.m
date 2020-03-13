function [rgbStates] = generateColorPanel(TotalNumber,SampleNumber,RGB)
    % the states of color for the sequence of tiles
    rgbStates = [];
    
    curState = RGB(1,:);
    rgbStates = cat(1,rgbStates,curState);
    
    
    for e = 2:(TotalNumber-1)
        
        if e <= SampleNumber + 1
            curState = RGB(2,:);
        end
        
        if e > SampleNumber 
            curState = RGB(3,:);
        end
        
        rgbStates = cat(1,rgbStates,curState);
    end
    curState = RGB(1,:);
    rgbStates = cat(1,rgbStates,curState);
end