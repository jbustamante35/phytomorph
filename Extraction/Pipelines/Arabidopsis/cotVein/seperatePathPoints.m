function [bPoints,sPoints] = seperatePathPoints(seedPoint,pathPoints)
    if size(pathPoints,1) == 0
        bPoints = seedPoint;
        sPoints = [];
    else
        ep = size(pathPoints,1);
        bpi = [1 ep];
        bPoints = pathPoints(bpi,:);
        pathPoints(bpi,:) = [];
        sPoints = pathPoints;
    end
end