function [] = bpCluster_func1(Adj,pointSets,threshold)
    cnt = 1;
    region = {};
    for s = 1:size(pointSets.bPoints,1)
        source = pointSets.bPoints(s,:);
        [path,pathcost] = searchAlongNetwork(Adj,pointSets.sPoints,pointSets.bPoints,source);
        gidx = find(pathcost < threshold);
        if numel(gidx) >= 2
            if (numel(region) < cnt);region{cnt} = [];end
            for e = 1:numel(gidx)
                region{cnt} = [region{cnt};pointSets.sPoints(path{gidx(e)},:)];
            end
        end
        cnt = cnt + 1;
    end
end