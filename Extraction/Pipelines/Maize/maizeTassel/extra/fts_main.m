%% sample the image and mask



fT = @(P,gD)transformFromGradient(gD,P);

% make sample domain
[x1,x2] = ndgrid(linspace(-domainW,domainW,domainW_num),linspace(-domainL,domainL,domainL_num));
szX = size(x1);
x = [x1(:) x2(:) ones(size(x1(:)))];

samI = [];
samM = [];

close all
tic
parfor e = 1:size(skel,1)
    tic
    samplePoint = skel(e,:);
    sampleP = funcP(samplePoint,fT);
    sampleP.getT(globG);
    
    samI(:,:,:,e) = sampleP.getF(I,x,szX);
    samM(:,:,e) = sampleP.getF(M,x,szX);
    toc
    
    %{
    imshow(samI(:,:,:,e),[]);
    title(num2str(e))
    drawnow
    %}
end
toc
%% make mask1
for e = 1:size(samM,3)
    tmpM = samM(:,:,e) > .5;
    hsz = hsize(tmpM);
    tmpM = imfill(~tmpM,hsz+.5,8) & tmpM;
    tmpR = regionprops(logical(tmpM),'MajorAxis','MinorAxis','Area','Perimeter','ConvexArea','Eccentricity');
    fm(e,:) = [tmpR.Area,tmpR.MajorAxisLength,tmpR.MinorAxisLength,tmpR.Eccentricity,tmpR.ConvexArea,tmpR.Perimeter];
    samM1(:,:,e) = tmpM;
    e
end
%% 
grps = 5;
zfm = zscore(fm);
options = statset('Display','iter');
gm = fitgmdist(zfm,grps,'Options',options);
kidx = gm.cluster(zfm);
close all
imshow(I,[]);hold on
for e = 1:grps
    plot(skel(kidx==e,1),skel(kidx==e,2),CL{e});
    ng(e) = sum(kidx==e);
end
%% features for images
for e = 1:size(samM,3)
    tmpM = samM(:,:,e);
    R = regionprops(logical(tmpM),'MajorAxis','MinorAxis','Area','Perimeter');
    f(e,1) = sum(tmpM(:));
end
%% page through masks
close all
for e = 1:size(samM,3)
    if f(e,1) < 1200
        out = flattenMaskOverlay(samI(:,:,:,e),logical(samM(:,:,e)));
        imshow(out,[]);
        drawnow
    end
end
%%
%%
domainW = 50;
domainL = 50;
domainW_num = domainW*2+1;
domainL_num = domainL*2+1;
% make sample domain
[x1,x2] = ndgrid(linspace(-domainW,domainW,domainW_num),linspace(-domainL,domainL,domainL_num));
szX = size(x1);
x = [x1(:) x2(:) ones(size(x1(:)))];

pointsToMatch = skel;
glf = @(X,M,m)((abs(X)>=1)*sign(X)*M) + ((abs(X)<1)*((M-m)*X+sign(X)*m));
sourcePoint = pointsToMatch(1000,:);
targetPoint = pointsToMatch(10000,:);

% convert RGB to gray and stack gradient
G = rgb2gray(I);
[g1,g2] = gradient(G);
globG = cat(3,G,g1,g2);
globG = cat(4,I,globG);

% make the function to create the frame at P
fT = @(P,gD)transformFromGradient(gD,P);

sourceP = funcP(sourcePoint,fT);
targetP = funcP(targetPoint,fT);

sourceF = sourceP.getF(globG,x,szX);
targetF = targetP.getF(globG,x,szX);

dX_range = [[50;0],[50;0]];
dR_range = [inf,0];
dS_range = [[2;.5],[2;.5]];

funcT = fwdT.makeTF(dX_range,dR_range,dS_range);
zp = dS_range(2,:).*-diff(dS_range,1,1).^-1;
pushT = funcT([0 0 .1 zp]);
sourceP.push(pushT)
sourceF2 = sourceP.getF(globG,x,szX);
close all
imshow([sourceF2,sourceF],[]);
close all
imshow([sourceF,targetF],[]);
%%

close all

tidx = randi(size(skel,1),1);
%%
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
%%
[sc,sv] = sort(c);
%%
close all
for e = 1:50
    %sidx = 81;
    sidx = sv(e);
    sourcePoint = pointsToMatch(sidx,:);
    sourceP = funcP(sourcePoint,fT);
    sourceP.getT(globG);
    sourceF = sourceP.getF(globG,x,szX);
    
    dT = funcT(tp(sidx,:));
    sourceP_copy = sourceP.push(dT);
    sourceF_copy = sourceP_copy.getF(globG,x,szX);
    imshow([sourceF sourceF_copy targetF],[]);
    drawnow
end



