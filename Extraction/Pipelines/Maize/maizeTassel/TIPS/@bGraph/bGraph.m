classdef bGraph < basisT
    % auto update the subgraphs upon modifying the main graph
    % the named sets need to be copied over
    % events would work for this
    
    properties
        % threshold for interconnects between nodes
        dtThresh;
        % internal graph object
        matG;
        % triangulation object for distances
        DT;
        % simple paths
        simplePaths;
        
        % affine embedding basis
        aBasis;
        
        % full simple paths
        fullSimplePaths;
        % complex junction paths
        complexPaths;
        % image of graph
        image;
        
        
        % named subsets
        namedSets;
        % sub-graps
        subGraphs;
        % name subgraphs
        namedsGraphs;
    end
    
    methods
    
        function [obj] = bGraph(graphPoints,dtThresh)
            if (size(graphPoints,2) == 2)
                graphPoints = [graphPoints,ones(size(graphPoints,1),1)];
            end
            if nargin == 1;dtThresh=2^.5;end
            obj = obj@basisT(graphPoints);
            obj.dtThresh = dtThresh;
            obj.updateGuts();
            obj.namedSets = struct('name',{},'position',{});
            obj.namedsGraphs = namedCollection('subGraphs','map');
        end
        
        function [nm] = lookUpName(this,queryPoints)
            if isa(queryPoints,'bGraph');queryPoints=queryPoints.E;end
            
            if size(queryPoints,2) ~= 1
                % loop over all the points
                for e = 1:size(queryPoints,1)
                    nm(e) = find(all(bsxfun(@eq,queryPoints(e,:),this.E),2));
                end
            else
                nm = queryPoints;
            end
        end
        
        function [X] = lookUpPosition(this,nm)
            if size(nm,2) == 1
                X = this.E(nm,:);
            else
                X = nm;
            end
        end
        
       
        
        
        
        function [newG] = neighbors(this,nm)
            nidx = this.matG.neighbors(nm);
            nHOOD = this.lookUpPosition(nidx);
            newG = bGraph(nHOOD);
        end
        
        function [X,idx] = nearest(this,queryPoint)
            for e = 1:size(queryPoint,1)
                %idx(e) = this.DT.nearestNeighbor(queryPoint(e,:));
                %idx(e) = dsearchn(this.E,this.DT,queryPoint(e,:));
                idx(e) = dsearchn(this.E,queryPoint(e,:));
            end
            X = this.lookUpPosition(idx);
        end
        
        function [bool] = contains(this,X)
            bool = isempty(setdiff(X,this.E,'rows'));
        end

        function [newG] = remove(this,Y)
            if isa(Y,'bGraph');Y = Y.E;end
            if nargout == 0
                this.E = setdiff(this.E,Y,'rows');
                this.removeFromSubset('all',Y);
                this.updateGuts();
            else
                newG = copy(this);
                newG.remove(Y);
            end
        end
        
        function [ret] = glue(X,Y)
            ret = [];
            if isa(X,'bGraph') && isa(Y,'bGraph')
                ret = bGraph([X.E;Y.E]);
            elseif isa(X,'bGraph') && ~isa(Y,'bGraph')
                X.E = [X.E;Y];
                X.updateGuts();
                ret = X;
            elseif ~isa(X,'bGraph') && isa(Y,'bGraph')
                Y.E = [X;Y.E];
                Y.updateGuts();
                ret = Y;
            end
        end
     
        function [newG] = intersect(this,A)
            d = intersect(this.E,A.E,'rows');
            newG = bGraph(d);
        end
        
        function [position,nidx] = getNhood(this,position)
            nm = this.lookUpName(position);
            nidx = this.matG.neighbors(nm);
            position = this.lookUpPosition(nidx);
        end
        
        function [distances,targets,graphR] = shortestpathtree(this,source,targets)
            if isa(source,'bGraph');source = source.E;end
            if isa(targets,'bGraph');targets = targets.E;end
            if isa(source,'char');source = this.namedSets(this.findSubset(source)).position;end
            if isa(targets,'char');targets = this.namedSets(this.findSubset(targets)).position;end
            graphR = [];distances = [];
            if ~isempty(targets)
                targets = this.lookUpName(targets);
                source = this.lookUpName(source);
                [graphR,distances] = this.matG.shortestpathtree(source,targets);
            end
        end
        
        function [pathX,pathIDX,B] = shortestpath(this,source,targets)
            try
                if isa(source,'bGraph');source = source.E;end
                if isa(targets,'bGraph');targets = targets.E;end
                if isa(source,'char');source = this.namedSets(this.findSubset(source)).position;end
                if isa(targets,'char');targets = this.namedSets(this.findSubset(targets)).position;end

                targets = this.lookUpName(targets);
                source = this.lookUpName(source);

                [pathIDX,B] = this.matG.shortestpath(source,targets);
                pathX = this.E(pathIDX,:);
            catch ME
                ME
            end
        end
        
        function [routes] = traceSimpleRoutes(this)
            routes = [];
            tmpG = this.namedsGraphs('simpleSegments');
            for e = 1:numel(tmpG)
                fprintf(['start:simple-analysis:' num2str(e) ':' num2str(numel(tmpG)) '\n']);
                sub1 = tmpG(e).getNamedSubset('routeEndPoints');
                sub2 = tmpG(e).getNamedSubset('branchPoints');
                if (size(sub1,1) == 2) && (size(sub2,1) == 0)
                    fprintf(['start:trace:' num2str(e) ':' num2str(numel(tmpG)) '\n']);
                    source = sub1(1,:);
                    target = sub1(2,:);
                    pathX = this.shortestpath(source,target);
                    routes = [routes,bCurve(pathX,this)];
                    fprintf(['end:trace:' num2str(e) ':' num2str(numel(tmpG)) '\n']);
                elseif (size(sub1,1) == 1) && (size(sub2,1) == 0)
                    pathX = sub1;
                    routes = [routes,bCurve(pathX,this)];
                end
                fprintf(['start:simple-analysis:' num2str(e) ':' num2str(numel(tmpG)) '\n']);
            end
            this.simplePaths = routes;
        end
        
        function [routes] = traceXtoY(this,X,Y)
            X = this.getNamedSubset(X);
            Y = this.getNamedSubset(Y);
            routes = [];
            for s = 1:size(X,1)
                fprintf(['start:simple-analysis:' num2str(s) ':' num2str(size(X,1)) '\n']);
                d = bsxfun(@minus,Y,X(s,:));
                d = find(sum(d.*d,2).^.5 < 30);
                [distances,targets] = this.shortestpathtree(X(s,:),Y(d,:));
                [~,midx] = min(distances);
                pathX = this.shortestpath(X(s,:),Y(d(midx),:));
                routes = [routes,bCurve(pathX,this)];
                fprintf(['start:simple-analysis:' num2str(s) ':' num2str(size(X,1)) '\n']);
            end
        end
        
        function [routes] = connectSimplePaths(this)
            routes = [];
            nHoodTH = 7;
            X = this.getNamedSubset('terminalPoints');
            for e = 1:numel(this.simplePaths)
                fprintf(['start:simple-analysis:' num2str(e) ':' num2str(numel(this.simplePaths)) '\n']);
             
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace start to terminal
                strPoint = this.simplePaths(e).strPoint';
                d = bsxfun(@minus,X,strPoint);
                d = find(sum(d.*d,2).^.5 < nHoodTH);
                
                tmpSet = X(d,:);
               
                [distances,targets] = this.shortestpathtree(strPoint,tmpSet);
                [~,midx] = min(distances);
                [pathXstr,pathIDX] = this.shortestpath(strPoint,tmpSet(midx,:));
                pathXstr = flip(pathXstr,1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace start to terminal
                stpPoint = this.simplePaths(e).stpPoint';
                d = bsxfun(@minus,X,stpPoint);
                d = find(sum(d.*d,2).^.5 < nHoodTH);
                
                tmpSet = X(d,:);
                tmpSet = setdiff(tmpSet,pathXstr(1,:),'rows');
                
                % for loops - if tmpSet if empty then allow for 
                % double use of solutions
                if isempty(tmpSet)
                    tmpSet = pathXstr(1,:);
                end
                
                
                [distances,targets] = this.shortestpathtree(stpPoint,tmpSet);
                [~,midx] = min(distances);
                pathXstp = this.shortestpath(stpPoint,tmpSet(midx,:));
                
                if size(this.simplePaths(e).E,1) ~= 1
                    pathX = [pathXstr;this.simplePaths(e).E(2:end-1,:);pathXstp];
                else
                    pathX = [pathXstr;pathXstp(2:end,:)];
                end
                
                
                routes = [routes,bCurve(pathX,this)];
                fprintf(['start:simple-analysis:' num2str(e) ':' num2str(numel(this.simplePaths)) '\n']);
            end
            this.fullSimplePaths = routes;
            
        end
        
        function [complexPaths] = connectComplexJunctions(this)
            complexPaths = [];
            jGraph = this.namedsGraphs('complexJunctions');
            for g = 1:numel(jGraph)
                branchPoints = jGraph(g).getNamedSubset('branchPoints');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for each branch point
                for p = 1:size(branchPoints,1)
                    % focus on this branch point
                    strPoint = branchPoints(p,:);
                    % get the nhood
                    [nHOOD,nidx] = this.getNhood(strPoint);
                    
                    % for each neighbor
                    for n = 1:numel(nidx)
                        % current neighbor point
                        curNpoint = this.E(nidx(n),:);
                        % is the curN a terminal point
                        if  this.subsetContains('terminalPoints',curNpoint)

                            % the path is short and full or terrors
                            pathX = [strPoint;curNpoint];
                            complexPaths = [complexPaths,bCurve(pathX,this)];

                        % if neighbor is NOT terminal point
                        else
                        
                            % make graph of ::skeleton-(focus + nhood(focus) - nX)
                            toRemove = setdiff([nHOOD;strPoint],curNpoint,'rows');
                            % remove points from graphs
                            potentialPathsGraph = jGraph(g).remove(toRemove);
                            % find the source point
                            sourceID = potentialPathsGraph.lookUpName(curNpoint);
                            % trace to terminal points
                            [D,TARGETS] = potentialPathsGraph.shortestpathtree(sourceID,'terminalPoints');
                            [D,midx] = min(D);
                            try
                                % if the target is not Inf away
                                if ~isempty(D) && ~isinf(D)
                                    % trace path to closest
                                    [pathX] = potentialPathsGraph.shortestpath(sourceID,TARGETS(midx));
                                    pathX = [strPoint;pathX];
                                    %pathX = bsxfun(@plus,pathX,deltaX);
                                    complexPaths = [complexPaths,bCurve(pathX,this)];
                                end
                            catch ME
                                here = 1;
                            end


                        end
                    end
                end
            end
            this.complexPaths = complexPaths;
        end
        
        
        function [subGraph] = makeSubGraph(this,X)
            subGraph = bGraph(X);
            subGraph.namedSets = this.intersectSubset('all',X);
        end
        %%%%%%%%%
        % subgraphs
        %%%%%%%%%
        %{
        function [] = addSubGraph(this,X)
            X = this.lookUpPosition(X);
            subGraph = bGraph(X);
            % copies over the named subsets from the main (this) graph
            subGraph.namedSets = this.intersectSubset('all',X);
            % stack the new graph onto the set of subgraphs
            this.subGraphs = [this.subGraphs,subGraph];
            %
        end
        %}
        
        function [] = addSubGraph(this,nameIndex,X)
            if (size(X,2) == 2)
                X = [X,ones(size(X,1),1)];
            end
            subGraph = bGraph(X);
            subGraph.namedSets = this.intersectSubset('all',X);
            if isempty(this.namedsGraphs(nameIndex))
                this.namedsGraphs(nameIndex) = subGraph;
            else
                %this.namedsGraphs(nameIndex).data(end+1) = subGraph;
                tmp = this.namedsGraphs(nameIndex);
                tmp(end+1) = subGraph;
                this.namedsGraphs(nameIndex) = tmp;
            end
            %{
            X = this.lookUpPosition(X);
            subGraph = bGraph(X);
            % copies over the named subsets from the main (this) graph
            subGraph.namedSets = this.intersectSubset('all',X);
            % stack the new graph onto the set of subgraphs
            this.subGraphs = [this.subGraphs,subGraph];
            %
            %}
        end
        
        
        %%%%%%%%%
        % named subsets
        % this is done using structs as the container for the named
        % subgraphs - we programmed a named Collection solution and 
        % are choosing to use it later for the named subsets
        %%%%%%%%%
        function [] = addNamedSet(this,name,X)
            if size(X,2) == 2;X = [X,ones(size(X,1),1)];end
            tmp.name = name;
            if ~isempty(X)
                X = this.lookUpPosition(X);
            end
            tmp.position = X;
            this.namedSets = [this.namedSets tmp];
        end
        
        function [result] = subsetContains(this,name,X)
            idx = this.findSubset(name);
            
            X = this.lookUpPosition(X);
            Y = this.namedSets(idx).position;
            
            result = intersect(Y,X,'rows');
            result = ~isempty(result);
        end
        
        function [result] = removeFromSubset(this,name,X)
            if strcmp(name,'all')
                name = {this.namedSets.name};
            end
            
            if ~iscell(name);name = {name};end
            
            for e = 1:numel(name)
                idx = this.findSubset(name{e});

                X = this.lookUpPosition(X);
                Y = this.namedSets(idx).position;

                tmp = setdiff(Y,X,'rows');
                if nargout == 0
                    this.namedSets(idx).position = tmp;
                else
                    result(e).name = this.namedSets(idx).name;
                    result(e).position = tmp;
                end
            end
        end
        
        function [result] = intersectSubset(this,name,X)
            if isa(X,'bGraph');X = X.E;end
            if strcmp(name,'all')
                name = {this.namedSets.name};
            end

            if ~iscell(name);name = {name};end

            for e = 1:numel(name)
                idx = this.findSubset(name{e});

                X = this.lookUpPosition(X);
                Y = this.namedSets(idx).position;

                tmp = intersect(Y,X,'rows');
                if nargout == 0
                    this.namedSets(idx).position = tmp;
                else
                    result(e).name = this.namedSets(idx).name;
                    result(e).position = tmp;
                end
            end
                
        end
        
        function [idx] = findSubset(this,name)
            idx = find(strcmp({this.namedSets.name},name));
        end
        
        function [sub] = getNamedSubset(this,name)
            dumpFlag = false;
            if ~iscell(name);name = {name};dumpFlag = true;end
            for e = 1:numel(name)
                idx = this.findSubset(name{e});
                sub{e} = this.namedSets(idx).position;
            end
            if dumpFlag;sub = sub{1};end
        end
        
        function [] = updateGuts(obj)
            %obj.DT = delaunayn(obj.E);
            %obj.DT = delaunayTriangulation(obj.E);
            
            distM = pdist2(obj.E,obj.E);
            distMSK = distM <= obj.dtThresh;
            ADJ = distM.*distMSK;
            obj.matG = graph(ADJ,'OmitSelfLoops');
            
            
        end
        
        function [h] = plot(this,cl,h)
            if nargin == 2;h=figure;imshow(this.image,[]);hold on;else;figure(h);end
            toProcess = [this.fullSimplePaths,this.complexPaths];
            toProcess.plot(cl);
        end
        
    end
    
end