%% list out the contents for the return folder- needed
FilePath = '/mnt/spaldingdata/nate/joeTest/';
FileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG','png','csv','txt'};
FileList = dig(FilePath,FileList,FileExt,1);
%% dig for files in pile - dig over local pile
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/tassel/UIUC_2019_Tassel_images/';
rawFileList = {};
FileExt = {'jpg'};
rawFileList = fdig(FilePath,FileList,FileExt,1);
%% reduce list to file name
orgName = {};
for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    [kvp] = findKVP(nm,'originalName');
    orgName{e} = getValue(kvp);
end
%% reduce original list
rawNM = {};
for e = 1:numel(rawFileList)
    [pth,nm,ext] = fileparts(rawFileList{e});
    rawNM{e} = nm;
end
%% match
clear dataCluster
for e = 1:numel(rawNM)
    idx = find(strcmp(orgName,rawNM{e}));
    tmp = FileList(idx);
    
    dataCluster(e).imageName = rawFileList{e};
    
    for i = 1:numel(tmp)
        [pth,nm,ext] = fileparts(tmp{i});
        if contains(nm,'maskNumber') && ...
                ~contains(nm,'processedImage') && ...
                ~contains(nm,'csvPhenoList')
            [kvp] = findKVP(nm,'maskNumber');
            maskIndex = str2num(getValue(kvp));
            dataCluster(e).maskSet{maskIndex} = tmp{i};
        elseif contains(nm,'processedImage') &&  ~contains(nm,'csvPhenoList')
            [kvp] = findKVP(nm,'maskNumber');
            maskIndex = str2num(getValue(kvp));
            dataCluster(e).processedSet{maskIndex} = tmp{i};
        elseif contains(nm,'csvPhenoList')
            [kvp] = findKVP(nm,'maskNumber');
            maskIndex = str2num(getValue(kvp));
            dataCluster(e).csvSet{maskIndex} = tmp{i};
        end
    end
end
%% load some data from csv for each dataCluster
for e = 1:numel(dataCluster)
    dataCluster(e).numCSV = numel(dataCluster(e).csvSet);
    cyverseList = {};
    errList = {};
    phenoVectors = [];
    for c = 1:numel(dataCluster(e).csvSet)
        try
            d = readtext(dataCluster(e).csvSet{c});
            cyverseList{c} = d{1};
            errList{c} = d{end};

            d(1) = [];
            d(end) = [];
            d = cell2mat(d);
            phenoVectors = [phenoVectors;d];
        catch ME
            
            cyverseList{c} = '';
            errList{c} = 'POST_PROCESS: Error reading';
            d = NaN*ones(1,8);
            phenoVectors = [phenoVectors;d];
        end
    end
    
    dataCluster(e).cyverseList = cyverseList;
    dataCluster(e).errList = errList;
    dataCluster(e).phenoVectors = phenoVectors';
    dataCluster(e).imageNameVec = e*ones(1,size(phenoVectors,1));
    dataCluster(e).maskNameVec = 1:size(phenoVectors,1);
    e
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load a mask for graph analysis
M = imread(dataCluster(5).maskSet{3});
%% start processing
smallAreaFilter = 70;
% fill small holes
M = bwareaopenSmallHoles(M,smallAreaFilter);
skel = bwmorph(M,'skel',inf);
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
% NOTE: FOR CHUNCK OF BRANCH POINTS - WE HAVE OR OPERATOR ABOVE
% BECAUSE INTERSECTION OF ENDPOINTS AND BRANCH POINTS IS NOT EMPTY
% SETE DIF IS NEEDED

%% goal: connect skeleton via branch points
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
junctionComplex = logical(skel - segs);
% find the endpoints of the segments
epSegs = bwmorph(segs,'endPoints');
[epSegsX(:,2),epSegsX(:,1)] = find(epSegs);
% create the graph
tasselGraph = bGraph(spX);
tasselGraph.image = M;
tasselGraph.addNamedSet('branchPoints',bpX);
tasselGraph.addNamedSet('endPoints',epX);
tasselGraph.addNamedSet('terminalPoints',[bpX;epX]);
tasselGraph.addNamedSet('routeEndPoints',epSegsX);
% generate the viewing image
RGB = cat(3,M,segs,junctionComplex);
%% add simple segments to the graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = regionprops(segs,'PixelIdxList','Image','BoundingBox');
for e = 1:numel(R)
    sub = [];
    [sub(:,2),sub(:,1)] = ind2sub(size(M),R(e).PixelIdxList);
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
    [sub(:,2),sub(:,1)] = ind2sub(size(M),eR(e).PixelIdxList);
    tasselGraph.addSubGraph('complexJunctions',sub);
end
%% trace and connect the simple routes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simpleRoutes = tasselGraph.traceSimpleRoutes();
fullSimpleRoutes = tasselGraph.connectSimplePaths();
complexPaths = tasselGraph.connectComplexJunctions();
%% remove end segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
imshow(double(RGB),[]);hold on
toProcess = fullSimpleRoutes;
CL = {'c','g','b'};
for l = 1:3
    nep = [];
    for e = 1:numel(toProcess)
        nep(e) = size(toProcess(e).getNamedSubset('endPoints'),1);
    end
    level{l} = toProcess(nep < 1);
    toProcess(nep<1) = [];
    level{l}.plot(CL{l});
    drawnow
end
simpleRoutes.plot('b');
complexPaths.plot('k');
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
imshow(double(RGB),[]);hold on
fullSimpleRoutes.plot('g');
simpleRoutes.plot('b');
complexPaths.plot('k');
[1603 3350];
[2050 4172];