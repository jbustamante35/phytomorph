classdef branchRegion < matlab.mixin.Heterogeneous & handle

    properties
    
        bPoints;
        sPoints;

    end

    methods

        function [obj] = branchRegion(bPoints,sPoints,cGraph,th)
            try
                obj.bPoints = bPoints;
                obj.sPoints = sPoints;
                cnt = 1;
                while cnt <= size(obj.bPoints,1)
                    seedPoint = obj.bPoints(cnt,:);
                    [path,pathcost] = cGraph.searchForBranchPoints(seedPoint);
                    fidx = find(pathcost < th);

                    for r = 1:numel(fidx)
                        pathPoints = cGraph.sPoints(path{fidx(r)},:);
                        [bPoints,sPoints] = seperatePathPoints(seedPoint,pathPoints);
                        if ~obj.inRegion(bPoints)
                            obj.appendPoints(bPoints,sPoints);
                        end
                    end
                    cnt = cnt + 1;
                end
            catch ME
                ME
            end
        end

        function [] = plot(obj)
            tPoints = [obj.bPoints;obj.sPoints];
            if size(tPoints,1) > 1
                K = convhull(tPoints(:,1)+rand(size(tPoints,1),1),tPoints(:,2)+rand(size(tPoints,1),1));
                plot(tPoints(K,2),tPoints(K,1),'r')
            else
                plot(tPoints(:,2),tPoints(:,1));
            end
        end

        function [b] = inRegion(obj,seedPoint)
            for s = 1:size(seedPoint,1)
                for e = 1:numel(obj)
                    b(e,s) = ~isempty(intersect([obj(e).bPoints;obj(e).sPoints],seedPoint(s,:),'rows'));
                end
            end
            b = all(b == 1,2);
        end

        function [] = appendPoints(obj,bPoints,sPoints)
            obj.bPoints = unique([obj.bPoints;bPoints],'rows','stable');
            obj.sPoints = unique([obj.sPoints;sPoints],'rows','stable');
        end


    end

end