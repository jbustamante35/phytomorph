function [skeleton,pointSets,Adj] = extractCotVeins(veinMask,gapClose,lengthFilter,smallHoleFilter,snipAmount)


    % fill small holes
    [veinMask] = fillSmallHoles(veinMask,smallHoleFilter);
    % create the skeleton
    skeleton = bwmorph(veinMask,'skeleton',inf);
    % find the branch points
    bp = bwmorph(skeleton,'branchpoints');
    % find the end points
    ep = bwmorph(skeleton,'endpoints');



    % find the branch,end and skeleton points
    [bPoints(:,1),bPoints(:,2)] = find(bp);
    [ePoints(:,1),ePoints(:,2)] = find(ep);
    [sPoints(:,1),sPoints(:,2)] = find(skeleton);
    % create the adj matrix
    Adj = Radjacency(sPoints',2^.5);



    % remove endpoints less then X to the nearest branch point
    pathcost = [];
    path = {};
    for stri = 1:size(ePoints,1)
        source = ePoints(stri,:);
        [tmpPath,tmpPathcost] = searchAlongNetwork(Adj,sPoints,bPoints,source);
        for stpi = 1:numel(tmpPath)
            path{stri,stpi} = tmpPath{stpi};
        end
        pathcost = [pathcost;tmpPathcost];
        fprintf(['Done with end point:' num2str(stri) ':' num2str(size(ePoints,1)) '\n'])
    end
    [pcost,pidx] = min(pathcost,[],2);
    ridx = find(pcost < snipAmount);

    % trim the 
    for r = 1:numel(ridx)
        rPath = path{ridx(r),pidx(ridx(r))};
        skeleton = removeBranch(skeleton,sPoints,rPath,bPoints);
    end
    


    
    skeleton = bwmorph(skeleton,'skeleton',inf);
    % remove objects with area less the length filter
    skeleton = bwareaopen(skeleton,lengthFilter);
    % re-skeletonize
    skeleton = bwmorph(skeleton,'skeleton',inf);
    % fill small holes
    skeleton = fillSmallHoles(skeleton,smallHoleFilter);
    % re-skeletonize
    skeleton = bwmorph(skeleton,'skeleton',inf);
    % close the skeleton
    skeleton = imclose(skeleton,strel('disk',gapClose,0));
    % re-skeletonize
    skeleton = bwmorph(skeleton,'skeleton',inf);


    pointSets.sPoints = sPoints;
    pointSets.ePoints = ePoints;
    pointSets.bPoints = bPoints;
    
end