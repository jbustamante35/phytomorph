function [mainVeinPath,mainVeinLength] = getMainVein(Adj,pointSets,bottomPoint,topPoint,skelPointsBool)
    mainVeinPath = [];
    mainVeinLength = 0;
    % find the main vein
    if skelPointsBool
        source = bottomPoint;
        target = topPoint;
        [mainVeinPath,mainVeinLength] = searchAlongNetwork(Adj,pointSets.sPoints,target,source);
        mainVeinPath = pointSets.sPoints(mainVeinPath{1},:);
    end
end