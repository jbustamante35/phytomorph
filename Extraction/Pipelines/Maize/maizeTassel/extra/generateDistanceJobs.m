function [] = generateDistanceJobs(sourceRGBName)

domainData.xData;
domainData.yData;



sourceData.rgb_fileName;
sourceData.mask_fileName;
sourceData.pixelList;


targetData.rgb_fileName;
targetData.mask_fileName;
targetData.pixelList;



M = double(M);
tic
targetPoint = pointsToMatch(tidx,:);
clear tp;

targetP = funcP(targetPoint,fT);
targetP.getT(globG);

targetF = targetP.getF(globG,x,szX);
targetM = targetP.getF(M,x,szX);
%targetM = targetM > .5;
[tx(:,1),tx(:,2)] = find(targetM);

dX_range = [[50;0],[50;0]];
dR_range = [inf,0];
dS_range = [[2;.5],[2;.5]];

funcT = fwdT.makeTF(dX_range,dR_range,dS_range);



for r = 1:100
    
    %sidx = randi(size(skel,1),1);
    sidx = r;
    
    sourcePoint = pointsToMatch(sidx,:);
    
    
    
    %targetPoint = sourcePoint + [10 10];

    
    sourceP = funcP(sourcePoint,fT);
    sourceP.getT(globG);


    %sourceF = sourceP.getF(globG,x,szX);
    %sourceM = sourceP.getF(M,x,szX);


    dT = [0 0 .1 zp];

    df = @(t)contrast1(t,funcT,sourceP,globG,M,targetF,targetM,x,szX);
    ops = optimset('Display','none');
    tic
    [tp(r,:),c(r)] = fminsearch(df,dT,ops);
    toc
    
    ops = optimoptions('particleswarm','Display','none');
    tic
    xsol = particleswarm(df,5,-1*ones(5,1),1*ones(5,1),ops);
    toc
end
toc
end