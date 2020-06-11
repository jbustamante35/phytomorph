function [newTasselGraph] = trimTerminalBranches(tasselGraph,lengthTH)
    if nargin == 1;lengthTH = inf;end
    toProcess = [tasselGraph.fullSimplePaths,tasselGraph.complexPaths];
    for e = 1:numel(toProcess)
        nep(e) = size(toProcess(e).getNamedSubset('endPoints'),1);
        len(e) = toProcess(e).length;
    end
    %rm = (nep < 1 & len < lengthTH);
    rm = (nep >= 1);
    toProcess(rm) = [];
    %% stack the points
    skelStack = [];
    for e = 1:numel(toProcess)
        skelStack = [skelStack;toProcess(e).E];
    end
    
    allBranchPoints = tasselGraph.getNamedSubset('branchPoints');
    %allEndPoints = tasselGraph.getNamedSubset('endPoints');
    %allBranchPoints = setdiff(allBranchPoints,allEndPoints,'rows');
    %skelStack = setdiff(skelStack,allBranchPoints,'rows');
    skelStack = [skelStack;allBranchPoints];
    
    skelIDX = sub2ind(size(tasselGraph(1).image),skelStack(:,2),skelStack(:,1));
    newSkel = zeros(size(tasselGraph(1).image));
    newSkel(skelIDX) = 1;
    newSkel = logical(newSkel);
    newSkel = bwmorph(newSkel,'skel',inf);
    newTasselGraph = makeTasselGraph(logical(newSkel));
end