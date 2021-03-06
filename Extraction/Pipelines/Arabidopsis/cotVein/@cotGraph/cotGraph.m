classdef cotGraph < handle
    
    properties
        skeletonI;
        sPoints;
        bPoints;
        ePoints;
        Adj;

        topPoint;
        botPoint;
        bRegions;

        branchRegions;

        mainVein;
        mainVeinLength;

        % non-main branch points
        nonMainBranchPoints;
        % main branch points
        mainBranchPoints;
        % main branch points - top abd bot
        midMainBranchPoints;
     
    end

    methods

        function [obj] = cotGraph(Adj,sPoints,ePoints,bPoints,skeletonI)
            obj.Adj = Adj;
            obj.sPoints = sPoints;
            obj.ePoints = ePoints;
            obj.bPoints = bPoints;
            obj.skeletonI = skeletonI;

            obj.botPoint;
            obj.topPoint;
            obj.mainVein;
        end

        function [] = decorateGraph(obj)
            obj.getTopPoint();
            obj.getBotPoint();
            obj.getMidVein();
            obj.divideBranchPoints();
        end

        function [] = getTopPoint(obj)
            branchPointsBool = ~isempty(obj.bPoints);
            edgePointsBool = ~isempty(obj.ePoints);
            skelPointsBool = ~isempty(obj.sPoints);
            % find the top and bottom point
            obj.topPoint = findTopBranchPoint(obj.Adj,obj,skelPointsBool,branchPointsBool,false);
        end

        function [] = getPointSets(obj)
            pointSets.sPoints = obj.sPoints;
        end

        function [] = getBotPoint(obj)
            branchPointsBool = ~isempty(obj.bPoints);
            edgePointsBool = ~isempty(obj.ePoints);
            skelPointsBool = ~isempty(obj.sPoints);
            obj.botPoint = findBottomEndPoint(obj.Adj,obj,skelPointsBool,edgePointsBool);
        end

        function [] = divideBranchPoints(obj)
            % non-main branch points
            obj.nonMainBranchPoints = setdiff(obj.bPoints,obj.mainVein,'rows');
            % main branch points
            obj.mainBranchPoints = setdiff(obj.bPoints,obj.nonMainBranchPoints,'rows');
            % main branch points - top abd bot
            obj.midMainBranchPoints = setdiff(obj.mainBranchPoints,obj.botPoint,'rows');
            obj.midMainBranchPoints = setdiff(obj.midMainBranchPoints,obj.topPoint,'rows');
        end

        function [] = getMidVein(obj)
            branchPointsBool = ~isempty(obj.bPoints);
            edgePointsBool = ~isempty(obj.ePoints);
            skelPointsBool = ~isempty(obj.sPoints);
            % trace the main vein
            [obj.mainVein,obj.mainVeinLength] = getMainVein(obj.Adj,obj,obj.botPoint,obj.topPoint,skelPointsBool);
        end

        function [] = generateRegions(obj,th)
            targetList = obj.mainBranchPoints;
            for e = 1:size(obj.mainBranchPoints,1)
                seedPoint = obj.mainBranchPoints(e,:);

                if ~isempty(obj.branchRegions)
                    if all(obj.branchRegions.inRegion(seedPoint)==0)
                        [path,pathcost] = searchAlongNetwork(obj.Adj,obj.sPoints,targetList,seedPoint);
                        fidx = find(pathcost < th);
                        for r = 1:numel(fidx)
                            pathPoints = obj.sPoints(path{fidx(r)},:);
                            [bPoints,sPoints] = seperatePathPoints(seedPoint,pathPoints);
                            if numel(obj.branchRegions) == 0
                                newRegion = branchRegion(bPoints,sPoints,obj,th);
                                obj.branchRegions = [obj.branchRegions;newRegion];
                            else
                                check = obj.branchRegions.inRegion(bPoints);
                                if all(check==0)
                                    newRegion = branchRegion(bPoints,sPoints,obj,th);
                                    obj.branchRegions = [obj.branchRegions;newRegion];
                                end
                            end
                        end
                    end
                else
                    [path,pathcost] = searchAlongNetwork(obj.Adj,obj.sPoints,targetList,seedPoint);
                    fidx = find(pathcost < th);
                    for r = 1:numel(fidx)
                        pathPoints = obj.sPoints(path{fidx(r)},:);
                        [bPoints,sPoints] = seperatePathPoints(seedPoint,pathPoints);
                       
                       
                        if numel(obj.branchRegions) == 0
                            newRegion = branchRegion(bPoints,sPoints,obj,th);
                            obj.branchRegions = [obj.branchRegions;newRegion];
                        else
                            check = obj.branchRegions.inRegion(bPoints);
                            if all(check==0)
                                newRegion = branchRegion(bPoints,sPoints,obj,th);
                                obj.branchRegions = [obj.branchRegions;newRegion];
                            end
                        end
                    end
                end
            end
        end

        function [path,pathcost] = searchForBranchPoints(obj,seedPoint)
            targetList = obj.mainBranchPoints;
            [path,pathcost] = searchAlongNetwork(obj.Adj,obj.sPoints,targetList,seedPoint);
        end

        function [] = initFromMask(veinMask,gapClose,lengthFilter,smallHoleFilter,snipAmount)

            [skeleton,pointSets] = extractCotVeins(veinMask,gapClose,lengthFilter,smallHoleFilter,snipAmount);

        end

        function [] = searchAlongNetwork()
            [path,pathcost] = searchAlongNetwork(Adj,networkList,targetList,source);
        end

        function [] = findInNetwork()
            [q1] = findInNetwork(networkList,source)
        end

        function [set] = getSet(obj,idx)
            set = obj.sPoints(idx,:);
        end

        function [] = measureAngles(obj)
            here = 1;
        end

        function [] = view()
        
        end

    end

end