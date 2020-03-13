%% mount max folder local
remoteBase = '/iplant/home/mbraud/UIUC_2019_Tassel_Images/';
localBase = '/home/nate/mntCyverse/tasselMax/';
mountCyVerse(remoteBase,localBase)
%sudo ./mountCyverse mbraud/UIUC_2019_Tassel_Images/ /home/nate/mntCyverse/tasselMax/
%% dig for image filess
FilePath = '/home/nate/mntCyverse/tasselMax/';
FileList = {};
FileExt = {'jpg'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
TIPS_noBackground(FileList{1},'./outputT/');
%%
oPath = '/mnt/spaldingdata/nate/tasselMasks/';
mmkdir(oPath);
parfor e = 1:numel(FileList)
    makeImageMask(FileList{e},oPath,10);
end
%%
samW = 50;
samN = samW*2 + 1;
% make sample square
[sam1,sam2] = ndgrid(linspace(-samW,samW,samN),linspace(-samW,samW,samN));
samSZ = size(sam1);
SAM = [sam1(:),sam2(:),ones(size(sam1(:)))];

% for each file
for e = 1:numel(FileList)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the file name
    [~,nm] = fileparts(FileList{e});
    % make the mask name
    maskName = [oPath nm '.tif'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    I = imread(FileList{e});
    I = double(I)/255;
    % read the mask
    M = imread(maskName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process the mask - get crop box
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make logical and get crop box
    M = logical(M);
    M = bwlarge(M);
    R = regionprops(M,'BoundingBox');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % crop the image and the mask
    I = imcrop(I,R(1).BoundingBox);
    M = imcrop(M,R(1).BoundingBox);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make overlay
    out = flattenMaskOverlay(I,M);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make skeleton, endpoints, and branch points list
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make skeleton
    skeleton = bwmorph(M,'skeleton',inf);
    % find skeleton points
    skel = [];
    [skel(:,2),skel(:,1)] = find(skeleton);
    
    % find end points
    ep = [];
    endpoints = bwmorph(skeleton,'endpoints');
    [ep(:,2),ep(:,1)] = find(endpoints);
    [epi] = findIndexInSpace(skel,ep);
    
    % find branch points
    bp = [];
    branchpoints = bwmorph(skeleton,'branchpoints');
    [bp(:,2),bp(:,1)] = find(branchpoints);
    [bpi] = findIndexInSpace(skel,bp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make distance transform
    DT = double(bwdist(~M));
    DT = imfilter(DT,fspecial('gaussian',[31 31],7),'replicate');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    [gDT1,gDT2] = gradient(DT);
    affine = zeros(1,3);
    affine(3) = 1;
    
    
    %M = cat(3,double(M),double(M),double(M));
    
    M = double(M);
    tmpI = zeros([samSZ size(I,3) size(skel,1)]);
    tmpM = zeros([samSZ size(M,3) size(skel,1)]);
   
    for p = 1:size(skel,1)
        n1 = gDT1(skel(p,2),skel(p,1));
        n2 = gDT2(skel(p,2),skel(p,1));
        t1 = -n2;
        t2 = n1;
        
        NOR = [n1 n2];
        NOR = NOR / norm(NOR);
        TAN = [t1 t2];
        TAN = TAN / norm(TAN);
        
        DIS = [skel(p,1) skel(p,2)];
        T = [[TAN',NOR',DIS'];affine];
        tSAM = mtimesx(T,SAM,'T')';
        
       
        tmp = ba_interp2(I,tSAM(:,1),tSAM(:,2));
        tmp = reshape(tmp,[samSZ size(I,3)]); 
        tmpI(:,:,:,p) = tmp;
        
        
       
        tmp = ba_interp2(M,tSAM(:,1),tSAM(:,2));
        tmp = reshape(tmp,[samSZ size(M,3)]);
        tmpM(:,:,:,p) = tmp;
        
        p
        size(skel,1)
    end
    
    
    %compute features
    NPP = 100;
    pi = zeros(size(tmpM,4),NPP);
    pSpace = linspace(0,norm(2*[samW samW]),NPP);
    parfor p = 1:size(tmpM,4)
        m = tmpM(:,:,:,p);
        i = tmpI(:,:,:,p);
        mn = m(:) / sum(m(:));
        mx = [];
        [mx(:,2),mx(:,1)] = find(m);
        Ex(p,:) = mtimesx(mn,'T',SAM);
        cp = Ex(p,:);
        cp = cp(1:2);
        m = bsxfun(@minus,mx,cp);
        dm = sum(m.*m,2).^.5;
        pi(p,:) = ksdensity(dm,pSpace);
        p
    end
    
    
    pi = bsxfun(@times,pi,sum(pi,2).^-1);
    en = -sum(log(pi+eps).*pi,2);
    [iS,iC,iU,iE,iL,iERR,iLAM] = PCA_FIT_FULL(pi,2);
    w = sweepPCA(iC,iE,iU,3*iLAM(1).^.5,1,5);
    [N] = hist3(iC,[1000 1000]);
    N = imfilter(N,fspecial('gaussian',[31 31],11),'replicate');
    plot(dX,iC(:,1),'.');
    
    
    
    dX = sum(Ex(:,1:2).*Ex(:,1:2),2).^.5;
    [sdX,sidx] = sort(dX);
    
    ng = 3;
    %tmpX = [dX iC(:,1)];
    tmpX = [dX en];
    tmpX = zscore(tmpX);
    kidx = kmeans(tmpX,ng);
    DC = [];
    CL = {'r.','b.','g.','k.','m.','c.','y.'};
    figure;
    for k = 1:ng
        fidx = find(kidx==k);
        subX = tmpX(fidx,:);
        plot(subX(:,1),subX(:,2),CL{k});
        hold on
        tmpU = mean(subX,1);
        
        MIN = min(subX,[],1);
        MAX = max(subX,[],1);
        tmpU = .5*[MAX - MIN] + MIN;
        
        subX = bsxfun(@minus,subX,tmpU);
        subX = sum(subX.*subX,2).^.5;
        [~,midx] = min(subX);
        
        DC = [DC , tmpI(:,:,:,fidx(midx))];
        
    end
    
    figure;
    imshow(DC,[]);
    hold on
    pause(1);
    drawnow
    xp = 10:size(tmpI,2):numel(CL)*size(tmpI,2);
    yp = 50*ones(size(xp));
    for p = 1:numel(CL)
        plot(xp(p),yp(p),CL{p});
    end
    
    % show clustered skeleton
    figure;
    imshow(I,[]);
    hold on
    for k = 1:ng
        fidx = find(kidx==k);
        plot(skel(fidx,1),skel(fidx,2),CL{k})
    end
    
    
    for p = 1:100:numel(sidx)
        imshow(tmpI(:,:,:,sidx(p)),[]);
        drawnow
    end
    
    
    close all
    for p = 1:100:size(tmpV,4)
        imshow(tmpV(:,:,:,p),[]);
        drawnow
    end
    
    
    %{
    % make adj matrix
    d = pdist(skel);
    d = squareform(d);
    dm = d <= 2^.5;
    d = double(dm).*d;
    
    % make graph from d
    g = graph(d);
    
  
    disp = true;
    if disp
        imshow(out,[]);
        hold on;
        drawnow
    end
    traceD = [];
    for i = 1:numel(epi)
        for j = 1:numel(bpi)
            [traceP{i,j},traceD(i,j)] = g.shortestpath(epi(i),bpi(j));
        end
        [~,sidx] = min(traceD(i,:));
        if disp
            
            plot(ep(i,1),ep(i,2),'g.');
            plot(bp(sidx,1),bp(sidx,2),'r.');
            tidx = traceP{i,sidx};
            %plot(skel(:,1),skel(:,2),'b.')
            plot(skel(tidx,1),skel(tidx,2),'k');
            drawnow
        end
    end
    %}
    
    
    
    
    
    
    out = flattenMaskOverlay(I,M);
    imshow(out,[]);
    hold on
    plot(skel(:,1),skel(:,2),'b.');
    hold off
    drawnow
end