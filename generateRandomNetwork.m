function [totalStack] = generateRandomNetwork()
    widthMax = 3;
    heightMax = 5;
    percentComposition = .5;
    H = randi(heightMax,1);
    
    unit = @(x)x+1;

    totalStack = [];
    for h = 1:H
        W = randi(widthMax,1);
        currentWidthStack = [];
        for w = 1:W


            layerMaxHeight = 3;
            [tmpStack] = generateSegment(percentComposition,layerMaxHeight,unit,h,w);


            currentWidthStack = [currentWidthStack,tmpStack];

        end
    
        layerName = ['Width Stack_' num2str(h)];
        currentWidthStack = layerStack(layerName,currentWidthStack);


        totalStack = [totalStack;currentWidthStack];
    end
    
    
    layerName = ['Total Stack'];
    totalStack = layerStack(layerName,totalStack);
    totalStack.view()

end


function [tmpStack] = generateSegment(percentComposition,layerMaxHeight,unit,h,w)
    per = rand(1);
    if per > percentComposition
        [tmpStack] = generateLayer(unit,h,w);    
    else
        [tmpStack] = generateLayerStack(layerMaxHeight,unit,h,w);    
    end

end



function [tmpStack] = generateLayer(unit,h,w)
    layerName = ['branchStack_' num2str(h) '_' num2str(w)];
    tmpStack = functionHandleLayer(layerName,unit);
end




function [tmpStack] = generateLayerStack(layerMaxHeight,unit,h,w)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate a layerStack
    H = randi(layerMaxHeight);
    tmpStack = [];
    for h = 1:H
        layerName = ['randLayer_' num2str(h) '_' num2str(w)];
        layer = functionHandleLayer(layerName,unit);
        tmpStack = [tmpStack;layer];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    layerName = ['branchStack_' num2str(w)];
    tmpStack = layerStack(layerName,tmpStack);
end