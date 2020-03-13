function [topPoint] = findTopBranchPoint(Adj,pointSets,skelPointsBool,branchPointsBool,toPlot)
    topPoint = [NaN,NaN];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the top skeleton xor branch point
    if skelPointsBool
        topThresh = 50;
        [~,topIdx] = min(pointSets.sPoints(:,1));
        if branchPointsBool
            source = pointSets.sPoints(topIdx,:);
            target = pointSets.bPoints;
            [tmpPath,tmpCost] = searchAlongNetwork(Adj,pointSets.sPoints,target,source);
            [minCost,minIdx] = min(tmpCost);
            if minCost < topThresh
                topIdx = minIdx;
                topPoint = pointSets.bPoints(topIdx,:);
                topIdx = find(all(bsxfun(@eq,pointSets.sPoints,topPoint) == 1,2));
            end
            if toPlot;plot(pointSets.sPoints(topIdx,2),pointSets.sPoints(topIdx,1),'ko','MarkerSize',14,'MarkerFaceColor','y','LineWidth',2);end
        end
        topPoint = pointSets.sPoints(topIdx,:);
    end
end