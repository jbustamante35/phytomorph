function [tasselGraph] = makeTasselGraph(skel)
    % our method for branch points - make cross kernel
    ker = ones(3);ker(1,1) = 0;ker(1,3) = 0;ker(3,1) = 0;ker(3,3) = 0;
    % find branch points via nhood count
    bp = imfilter(double(skel),ker,'replicate') .* skel >= 4;
    % find end points
    ep = bwmorph(skel,'endPoints');
    % find branch points via binary morphology
    bp2 = bwmorph(skel,'branchPoints');
    % OR the branch points together
    bp = bp | bp2;
    % end-points, branch-points,skeleton-points
    epX = [];bpX = [];spX = [];
    epSegsX = [];
    ep = ep & ~ bp;
    % find the X locations
    [epX(:,2),epX(:,1)] = find(ep);
    [bpX(:,2),bpX(:,1)] = find(bp);
    [spX(:,2),spX(:,1)] = find(skel);
    %% goal: connect skeleton via branch points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define: branch point above
    % bp-nHood: the area near the branch point
    %           in this case is 8-connected
    % paths: list of (x,y)-pairs which include the branch point
    % and the points inbetween.
    % Solution: break into simple and complex paths.
    % Simple: paths are below as remove the nhood (dialated bp-map) from the skeleton
    % find the end points of the remaining segments and trace them.
    % Complex: create the nhood for each branch point.  If an nhood group
    % contains more than 1 branch point, then trace as follow:
    nhood = logical(imdilate(bp,ones(3)));
    % make segments image
    segs = (skel==1) & (nhood == 0);
    % find the endpoints of the segments
    epSegs = bwmorph(segs,'endPoints');
    [epSegsX(:,2),epSegsX(:,1)] = find(epSegs);
    % create the graph
    tasselGraph = bGraph(spX);
    tasselGraph.image = skel;
    tasselGraph.addNamedSet('branchPoints',bpX);
    tasselGraph.addNamedSet('endPoints',epX);
    tasselGraph.addNamedSet('terminalPoints',[bpX;epX]);
    tasselGraph.addNamedSet('routeEndPoints',epSegsX);
    %% add simple segments to the graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = regionprops(segs,'PixelIdxList','Image','BoundingBox');
    for e = 1:numel(R)
        sub = [];
        [sub(:,2),sub(:,1)] = ind2sub(size(skel),R(e).PixelIdxList);
        tasselGraph.addSubGraph('simpleSegments',sub);
    end
    %% add complex junctions to the graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    junctionComplex = logical(skel - segs);
    % find the complex junctions by counting the number of branch points in the
    eR = regionprops(junctionComplex,'PixelIdxList','Area','Image','Centroid','BoundingBox');
    numBranchPoints = [];
    for e = 1:numel(eR)
        numBranchPoints(e) = sum(bp(eR(e).PixelIdxList));
    end
    complexIDX = find(numBranchPoints >= 2);
    eR = eR(complexIDX);
    % add the complex junctions
    for e = 1:numel(eR)
        sub = [];
        [sub(:,2),sub(:,1)] = ind2sub(size(skel),eR(e).PixelIdxList);
        tasselGraph.addSubGraph('complexJunctions',sub);
    end
    %% trace and connect the simple routes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simpleRoutes = tasselGraph.traceSimpleRoutes();
    fullSimpleRoutes = tasselGraph.connectSimplePaths();
    complexPaths = tasselGraph.connectComplexJunctions();
end