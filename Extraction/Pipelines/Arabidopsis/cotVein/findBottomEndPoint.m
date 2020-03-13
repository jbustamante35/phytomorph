function [botPoint] = findBottomEndPoint(Adj,pointSets,skelPointsBool,edgePointsBool)
botPoint = [NaN,NaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if skelPointsBool
    % find the bottom endpoint
    if ~edgePointsBool
        [~,botIdx] = max(pointSets.sPoints(:,1));
        source = pointSets.sPoints(botIdx,:);
    else   
        [~,botIdx] = max(pointSets.ePoints(:,1));
        source = pointSets.ePoints(botIdx,:);
        [botIdx] = findInNetwork(pointSets.sPoints,source);
        source = pointSets.sPoints(botIdx,:);
    end
    %{
    % find the stem if can
    if branchPointsBool
        [tmpPath,tmpPathcost] = searchAlongNetwork(Adj,pointSets.sPoints,pointSets.bPoints,source);
        [~,midx] = min(tmpPathcost);
        maPath = pointSets.sPoints(tmpPath{midx},:);
        plot(maPath(:,2),maPath(:,1),'r.','MarkerSize',3);
    end
    %}
    botPoint = pointSets.sPoints(botIdx,:);
    %plot(pointSets.sPoints(botIdx,2),pointSets.sPoints(botIdx,1),'ko','MarkerSize',14,'MarkerFaceColor','k','LineWidth',2);
end