function [d] = myPointSetContrast(moving,fixed,d,trans,fixedPoints)

    edgeMoving = transformImage(moving,fixed,d,trans);
    [edgePointsM(:,1),edgePointsM(:,2)] = find(edgeMoving);

    d = pdist2(edgePointsM,fixedPoints);
    d1 = mean(min(d,[],1));
    d2 = mean(min(d,[],2));
    d = .5*(d1+d2);
end