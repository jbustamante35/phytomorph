FilePath = '/mnt/snapper/nate/Madison 19 research/';
FileList = {};
FileExt = {'tif'};
verbose = 1;
FileList = sdig(FilePath,FileList,FileExt,verbose);
%%
close all
dataStack = [];
dataStack_bin = [];
dataStack_dT = [];
samCNT = 1;
DMSZ = 20;
CH = 50000;
perCH = 50;
%CH = 100;
[n1,n2] = ndgrid(-DMSZ:DMSZ,-DMSZ:DMSZ);
sDomain = [n1(:) n2(:) ones(size(n1(:)))];
dataStack = zeros([size(n1) CH]);
dataStack_bin = zeros([size(n1) CH]);
dataStack_dT = zeros([size(n1) CH]);
dataStack_skel = zeros([size(n1) CH]);
dataStack_nor1 = zeros([size(n1) CH]);
dataStack_nor2 = zeros([size(n1) CH]);
s = 1;
%for s = 1:numel(FileList)
stopFlag = false;
disp = 0;
while ~stopFlag
    try
        
        fprintf(['Loading images....START\n']);
        %%%%%%%%%%%%%%%%%%%%%
        % load the images
        out = [];
        I = [];
        for e = 1:numel(FileList{s})
            I(:,:,e) = imread(FileList{s}{e});
        end
        I = uint8(I);
        % load the images
        %%%%%%%%%%%%%%%%%%%%
        fprintf(['Loading images....END\n']);
        %I = double(I)/255;

        %{
        % use monomodal
        [optimizer,metric] = imregconfig('monomodal');
        % multimodel does not work well
        [optimizer,metric] = imregconfig('multimodal');
        % align 1 to 2
        movingRegistered1 = imregister(I(:,:,2), I(:,:,1), 'rigid', optimizer, metric);
        % align 1 to 3
        movingRegistered2 = imregister(I(:,:,3), I(:,:,1), 'rigid', optimizer, metric);
        rI(:,:,1) = I(:,:,1);
        rI(:,:,2) = movingRegistered1;
        rI(:,:,3) = movingRegistered2;
        %}
        
        %%%%%%%%%%%%%%%%%%%%%
        % thresholds for binary object removal
        % area threshold
        areaThreshold = 20;
        eccentricityThreshold = .93;
        motionAreaThreshold = 10;
        dilateValue = 5;
        spurValue = 5;
        
        fprintf(['Making masks....START\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the back ground image
        mskF = [];
        for t = 1:size(I,3)
            tic
            fprintf(['(' num2str(t) ')-->start\n']);
            %%%%%%%%%%%%%%%%%%%%%
            % find the foreground
            tmp = double(I(:,:,t));
            tmp = imfilter(tmp,fspecial('average',61),'replicate');
            tmp = imdilate(tmp,strel('disk',71));
            tmp = imfilter(tmp,fspecial('average',61),'replicate');
            foreGround(:,:,t) = bindVec(double(I(:,:,t)) - tmp);
            msk1 = double(I(:,:,t)/255) < graythresh(double(I(:,:,t)/255));
            msk2 = double(foreGround(:,:,t)) < graythresh(foreGround(:,:,t));
            
            
            %%%%%%%%%%%%%%%%%%%%%
            % fill holes
            
            msk2 = imclose(msk2,strel('disk',dilateValue,0));
            msk2 = imfill(msk2,'holes');
            
            %%%%%%%%%%%%%%%%%%%%%
            % clear sides
            bottom = msk2(end,:);
            msk2(end,:) = 0;
            msk2 = imclearborder(msk2);
            msk2(end,:) = bottom;
            
            %%%%%%%%%%%%%%%%%%%%%
            % apply area threshold
            msk2 = bwareaopen(msk2,areaThreshold);
            
            %%%%%%%%%%%%%%%%%%%%%
            % clear sides
            
            out = flattenMaskOverlay(double(I(:,:,t))/255,msk1,.2,'r');
            out = flattenMaskOverlay(out,msk2,.2,'b');
            R = regionprops(msk2,'Perimeter','MajorAxisLength','MinorAxisLength',...
                                'Centroid','Area','Eccentricity','PixelIdxList');
            
            %%%%%%%%%%%%%%%%%%%%%
            % apply area threshold
            kidx1 = find([R.Eccentricity] > eccentricityThreshold);
            mskT = zeros(size(msk2));
            for obj = 1:numel(kidx1)
                mskT(R(kidx1(obj)).PixelIdxList) = 1;
            end
            
            
            mskF(:,:,t) = mskT;
            fprintf(['(' num2str(t) ')-->end @' num2str(toc) '\n']);
        end
        fprintf(['Making masks....END\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        fprintf(['Align Images....START\n']);
        rI = [];
        rI(:,:,1) = I(:,:,1);
        mskR(:,:,1) = mskF(:,:,1);
        [rI(:,:,2),mskR(:,:,2)] = steveAlign(I(:,:,[1 2]),mskF(:,:,2));
        [rI(:,:,3),mskR(:,:,3)] = steveAlign(I(:,:,[1 3]),mskF(:,:,3));
        fprintf(['Align Images....END\n']);
        
        wholeMask = mskR;
        
        baseDomain = regionprops(logical(any(wholeMask,3)),'Perimeter','MajorAxisLength','MinorAxisLength',...
                                'Centroid','Area','Eccentricity','PixelIdxList');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pixel change in intensity after reg
        dI1 = double(rI(:,:,2)) - double(rI(:,:,1));
        dI2 = double(rI(:,:,3)) - double(rI(:,:,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % analysis for 1 -> 2
        fidx = find(dI1 < 0);
        vec = bindVec(dI1(fidx));
        level = graythresh(vec);
        msk = find(vec < level);
        Motion_MASK1 = zeros(size(dI1));
        Motion_MASK1(fidx(msk)) = 1;
        Motion_MASK1 = imclearborder(Motion_MASK1);
        Motion_MASK1 = Motion_MASK1 .* wholeMask(:,:,2);
        Motion_MASK1 = bwareaopen(Motion_MASK1,motionAreaThreshold);
        
        moveDomain1 = regionprops(logical(Motion_MASK1),'Perimeter','MajorAxisLength','MinorAxisLength',...
                                'Centroid','Area','Eccentricity','PixelIdxList');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % analysis for 2 -> 3
        fidx = find(dI2 < 0);
        vec = bindVec(dI2(fidx));
        level = graythresh(vec);
        msk = find(vec < level);
        Motion_MASK2 = zeros(size(dI2));
        Motion_MASK2(fidx(msk)) = 1;
        Motion_MASK2 = imclearborder(Motion_MASK2);
        Motion_MASK2 = Motion_MASK2 .* wholeMask(:,:,3);
        Motion_MASK2 = bwareaopen(Motion_MASK2,motionAreaThreshold);
        moveDomain2 = regionprops(logical(Motion_MASK2),'Perimeter','MajorAxisLength','MinorAxisLength',...
                                'Centroid','Area','Eccentricity','PixelIdxList');
        
                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over objects in move1 domain to intersect with baseDomain
        kp = [];
        for baseDomainN = 1:numel(baseDomain)
            tmpKeep = [];
            X = baseDomain(baseDomainN).PixelIdxList;
            for moveDomainN = 1:numel(moveDomain1)
                Y = moveDomain1(moveDomainN).PixelIdxList;
                tmpKeep(moveDomainN) = ~isempty(intersect(X,Y));
            end
            kp(baseDomainN) = any(tmpKeep);
        end
        baseDomain = baseDomain(find(kp==1));
        % make a base domain image
        baseImage = zeros(size(rI,1),size(rI,2));
        for e = 1:numel(baseDomain)
            baseImage(baseDomain(e).PixelIdxList) = 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        Motion_MASK1 = Motion_MASK1 .* baseImage;
        Motion_MASK2 = Motion_MASK2 .* baseImage;
        wholeMask = bsxfun(@times,wholeMask,baseImage);
                            
        out = [];
        out(:,:,:,1) = cat(3,rI(:,:,1),rI(:,:,1),rI(:,:,1))/255;
        out(:,:,:,2) = flattenMaskOverlay(rI(:,:,2)/255,logical(Motion_MASK1));
        out(:,:,:,3) = flattenMaskOverlay(rI(:,:,3)/255,logical(Motion_MASK2));
        
        
        for t = 1:size(wholeMask,3)
            skel(:,:,t) = bwmorph(logical(wholeMask(:,:,t)),'skeleton',inf);
            skel(:,:,t) = bwmorph(skel(:,:,t),'spur',spurValue);
        end
        
        
        
        
        for t = 1:size(rI,3)
            
            out(:,:,:,t) = flattenMaskOverlay(out(:,:,:,t),logical(wholeMask(:,:,t)),.3,'b');
            mR{t} = regionprops(logical(wholeMask(:,:,t)),'Centroid','MajorAxisLength','Orientation','Area','PixelIdxList','BoundingBox');
            
            
            
            
            for obj = 1:numel(mR{t})
                objFilter = zeros(size(rI,1),size(rI,2));
                objFilter(mR{t}(obj).PixelIdxList) = 1;
                tmpS = skel(:,:,t) .* objFilter;
                mR{t}(obj).skelLength = sum(tmpS(:));
                
                
                [sy,sx] = find(tmpS);
                mR{t}(obj).skelPoints = [sx sy];
            end
            
        end
        
        fprintf(['Sampling Images...START'])
        for t = 1:size(out,4)
            
            fprintf(['(' num2str(t) ')-->start\n']);
             
             
            imshow(out(:,:,:,t),[]);
            hold on
            
            
            for obj = 1:numel(mR{t})
                fprintf(['\t(' num2str(obj) ')-->obj_start\n']);
                
                
                
                objFilter = zeros(size(rI,1),size(rI,2));
                objFilter(mR{t}(obj).PixelIdxList) = 1;
                objDomain = objFilter;
                [d2 d1] = find(objFilter);
                objFilter = bwmorph(objFilter,'skel',inf);
                objSkel = objFilter;
                dT = double(bwdist(objFilter));
                dT = imfilter(dT,fspecial('gaussian',[21 21],4),'replicate');
                [dX,dY] = gradient(dT);
                
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% sample the domain points
                %{
                ptToSample = [d1 d2];
                toSample = size(ptToSample,1);
                
                for p = 1:toSample
                    
                    
                    m = ptToSample(p,:);
                    tx = ba_interp2(dX,ptToSample(p,1),ptToSample(p,2));
                    ty = ba_interp2(dY,ptToSample(p,1),ptToSample(p,2));
                    tan = [tx ty];
                    tan = tan / norm(tan);
                    nor = [tan(2) -tan(1)];
                    aff = [[tan';0] [nor';0] [m';1]];
                    tmpDomain = (aff*sDomain')';
                    
                    sam = ba_interp2(rI(:,:,t),tmpDomain(:,1),tmpDomain(:,2));
                    sam = reshape(sam,size(n1));
                    c(:,p) = PCA_REPROJ_T(sam(:),mE,mU);
                    p
                    toSample
                    
                end
                %}
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% sample the skeleton points
              
                if samCNT <= CH
                    ptToSample = mR{t}(obj).skelPoints;
                    ptToSample = ptToSample(randperm(size(ptToSample,1)),:);
                    toSample = min(perCH,size(ptToSample,1));

                    for p = 1:toSample
                

                        m = ptToSample(p,:);
                        nx = ba_interp2(dX,ptToSample(p,1),ptToSample(p,2));
                        ny = ba_interp2(dY,ptToSample(p,1),ptToSample(p,2));
                        nor = [nx ny];
                        nor = nor / norm(nor);
                        tan = [nor(2) -nor(1)];
                        aff = [[tan';0] [nor';0] [m';1]];
                        iaff = inv(aff);
                        tmpDomain = (aff*sDomain')';

                        sam = ba_interp2(rI(:,:,t)/255,tmpDomain(:,1),tmpDomain(:,2));
                        sam = reshape(sam,size(n1));
                        %dataStack = cat(3,dataStack,sam);
                        dataStack(:,:,samCNT) = sam;
                        
                        % set view sample
                        samV = sam;
                        
                        sam = ba_interp2(objDomain,tmpDomain(:,1),tmpDomain(:,2));
                        sam = reshape(sam,size(n1));
                        %sam = imclose(sam,strel('disk',5,0));
                        %sam = imfill(sam,'holes');
                        %dataStack_bin = cat(3,dataStack_bin,sam);
                        dataStack_bin(:,:,samCNT) = sam;

                        sam = ba_interp2(dT,tmpDomain(:,1),tmpDomain(:,2));
                        sam = reshape(sam,size(n1));
                        %dataStack_dT = cat(3,dataStack_dT,sam);
                        dataStack_dT(:,:,samCNT) = sam;
                        
                        sam = ba_interp2(double(objSkel),tmpDomain(:,1),tmpDomain(:,2));
                        sam = reshape(sam,size(n1));
                        %dataStack_dT = cat(3,dataStack_dT,sam);
                        dataStack_skel(:,:,samCNT) = sam;
                        
                        sam1 = ba_interp2(dX,tmpDomain(:,1),tmpDomain(:,2));
                        sam1 = reshape(sam1,size(n1));
                        sam2 = ba_interp2(dY,tmpDomain(:,1),tmpDomain(:,2));
                        sam2 = reshape(sam2,size(n1));
                        samD = [sam1(:) sam2(:) zeros(size(sam1(:)))]';
                        samD = (iaff*samD)';
                        sam1 = reshape(samD(:,1),size(n1));
                        sam2 = reshape(samD(:,2),size(n1));
                        
                        dataStack_nor1(:,:,samCNT) = sam1;
                        dataStack_nor2(:,:,samCNT) = sam2;

                        samCNT = samCNT + 1;
                        
                        if samCNT < disp
                        
                            imshow(out(:,:,:,t),[]);
                            hold on
                            plot(tmpDomain(:,1),tmpDomain(:,2),'.')
                            plot(m(1),m(2),'r*')
                            imshow(samV,[]);
                            hold off
                            drawnow
                        end
                        
                        
                        fprintf('.');
                        if mod(p,30) == 0;fprintf('.\n');end
                    end
                else
                    stopFlag = true;
                    break
                end
                fprintf(['\n']);
                fprintf(['\t(' num2str(obj) ')-->obj_end\n']);
            end
            
            
            
            
            
            for obj = 1:numel(mR{t})
                
                
                
                ang = -mR{t}(obj).Orientation*pi/180;
                vec = [cos(ang) sin(ang)];
                vec1 = vec * mR{t}(obj).MajorAxisLength/2;
                vec2 = vec * mR{t}(obj).skelLength/2;
                
                
                plot(mR{t}(obj).skelPoints(:,1),mR{t}(obj).skelPoints(:,2),'g.','MarkerSize',2)
                
                
                m = mR{t}(obj).Centroid;
                quiver(m(1),m(2),vec1(1),vec1(2),'Color','r','LineWidth',1)
                quiver(m(1),m(2),-vec1(1),-vec1(2),'Color','r','LineWidth',1)
                
                quiver(m(1),m(2),vec2(1),vec2(2),'Color','b','LineWidth',1)
                quiver(m(1),m(2),-vec2(1),-vec2(2),'Color','b','LineWidth',1)
                plot(m(1),m(2),'b')
            end
            
           
            
            drawnow
            hold off
            pause(.2)
            %waitforbuttonpress
        end
        
        
    catch ME
        getReport(ME)
        
    end
    
    s = s + 1;
end
%% model all the data raw all
numC = 170;
dataStack_T = cat(4,dataStack,dataStack_bin,dataStack_dT,dataStack_skel,dataStack_nor1,dataStack_nor2);
dataStack_T = permute(dataStack_T,[1 2 4 3]);
szT = size(dataStack_T);
dataStack_T = reshape(dataStack_T,[prod(szT(1:3)) szT(4)]);
rm = any(isinf(dataStack_T) | isnan(dataStack_T),1);
dataStack_T(:,rm) = [];
szT(4) = size(dataStack_T,2);

[mU,mE,mL] = PCA_FIT_FULL_Tws(dataStack_T,numC);
mC = PCA_REPROJ_T(dataStack_T,mE,mU);
mS = PCA_BKPROJ_T(mC,mE,mU);

mS = reshape(mS,szT);
mS = permute(mS,[1 2 4 3]);


dataStack_T = reshape(dataStack_T,szT);
dataStack_T = permute(dataStack_T,[1 2 4 3]);



%% watch bin images
close all
for e = 1:30%size(dataStack_T,3)
    rawImage = [];
    simImage = [];
    for s = 1:size(dataStack_T,4)
        rawImage = cat(2,rawImage,bindVec(dataStack_T(:,:,e,s)));
        simImage = cat(2,simImage,bindVec(mS(:,:,e,s)));
    end
    imshow([rawImage;simImage],[]);
    drawnow
end
%% measure the raw data patches
close all
binTH = .1;
disp = 20;
NP = [];
MS = [];
CM = [];
EP_matrix = [];
domainArea = [];
frameImage = zeros(2*DMSZ+1,2*DMSZ+1);
for e = 1:4
    frameImage(:,1) = 1;
    frameImage = imrotate(frameImage,90); 
end

for e = 1:size(dataStack_T,3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skeleton measurements
    tmp = dataStack_T(:,:,e,4);
    tmp = logical(tmp > binTH);
    tmp = bwmorph(tmp,'skel',inf);
    tmp = bwmorph(tmp,'spur',3);
    ep = bwmorph(tmp,'endpoints');
    [e2,e1] = find(ep);
    eidx = find(ep);
    outOver1 = tmp;
    EP = [e1 e2];
    kpidx = frameImage(eidx) == 1;
    EP = EP(kpidx,:);
    [~,srtIDX] = sort(EP(:,2),'descend');
    EP = EP(srtIDX,:);
    NP(e) = size(EP,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % domain measurements
    tmp = logical(dataStack_T(:,:,e,2) > binTH);
    [d2,d1] = find(tmp);
    CM(e,:) = [mean(d1) mean(d2)];
    CM(e,:) = [21 21];
    domainArea(e) = sum(tmp(:));
    tmp1 = dataStack_T(:,:,e,end);
    tmp2 = dataStack_T(:,:,e,end-1);
    qX = [];
    qX(:,1) = ba_interp2(tmp1,EP(:,1),EP(:,2));
    qX(:,2) = ba_interp2(tmp2,EP(:,1),EP(:,2));
    qX_NOR = sum(qX.*qX,2).^-.5;
    qX = bsxfun(@times,qX,qX_NOR);
    
    qM = [];
    qM(:,1) = ba_interp2(tmp1,CM(e,1),CM(e,2));
    qM(:,2) = ba_interp2(tmp2,CM(e,1),CM(e,2));
    qM_NOR = sum(qM.*qM,2).^-.5;
    qM = bsxfun(@times,qM,qM_NOR);
    
    
    
    
    
    
    if NP(e) == 2
        vec1 = CM(e,:) - EP(1,:);
        vec2 = -(CM(e,:) - EP(2,:));
        vec1 = vec1 / norm(vec1);
        vec2 = vec2 / norm(vec2);
        MS(e) = vec1*vec2';
        EP_matrix = cat(3,EP_matrix,EP);
    elseif NP(e) == 1
        vec1 = CM(e,:) - EP(1,:);
        MS(e) = 1;
        EP_matrix = cat(3,EP_matrix,[EP;zeros(1,2)]);
    else 
        MS(e) = 0;
        EP_matrix = cat(3,EP_matrix,zeros(2));
    end
    
    
    if e < disp
        out = flattenMaskOverlay(dataStack_T(:,:,e,1),outOver1,1,'r');
        imshow(out,[]);
        hold on
        quiver(EP(:,1),EP(:,2),qX(:,1),qX(:,2),'m');
        quiver(EP(:,1),EP(:,2),-qX(:,1),-qX(:,2),'m');
        quiver(CM(e,1),CM(e,2),qM(1),qM(2),10,'m');
        quiver(CM(e,1),CM(e,2),-qM(1),-qM(2),10,'m');
          plot(EP(:,1),EP(:,2),'g.');
        plot(CM(e,1),CM(e,2),'b.');
        hold off
        drawnow
        axis([-5 50 -5 50])
        %waitforbuttonpress
    end
    %waitforbuttonpress
end
%% sub section the data
sidx = find(NP == 2);
%sidx = find(NP == 1);
numC = 170;
dataStack_T = cat(4,dataStack,dataStack_bin,dataStack_dT,dataStack_skel,dataStack_nor1,dataStack_nor2);
dataStack_T = permute(dataStack_T,[1 2 4 3]);
szT = size(dataStack_T);
dataStack_T = reshape(dataStack_T,[prod(szT(1:3)) szT(4)]);
rm = any(isinf(dataStack_T) | isnan(dataStack_T),1);
fprintf(['Removed TR:' num2str(sum(rm)) '\n']);
dataStack_T(:,rm) = [];

auxData1 = CM(sidx,:)';
aSZ1 = size(auxData1);

auxData2 = EP_matrix(:,:,sidx);
aSZ2 = size(auxData2);
auxData2 = reshape(auxData2,[prod(aSZ2(1:2)) aSZ2(3)]);


% perform sub section
dataStack_T = dataStack_T(:,sidx);
dataStack_T = [dataStack_T;auxData1;auxData2];


szT(4) = size(dataStack_T,2);

[mU,mE,mL] = PCA_FIT_FULL_Tws(dataStack_T,numC);
mC = PCA_REPROJ_T(dataStack_T,mE,mU);
mS = PCA_BKPROJ_T(mC,mE,mU);

% peal off auxData and pop from stack
mA1 = mS((end-5):(end-4),:);
mA1 = mA1';

mA2 = mS((end-3):end,:);
mA2 = reshape(mA2,aSZ2);
mS((end-5):end,:) = [];

% peal off data stack
dataStack_T((end-5):end,:) = [];

mS = reshape(mS,szT);
mS = permute(mS,[1 2 4 3]);
UQ = unique(dataLabel);


dataStack_T = reshape(dataStack_T,szT);
dataStack_T = permute(dataStack_T,[1 2 4 3]);

%% watch bin images - SUB SECTION
close all

subMS = MS(sidx);
[~,ssidx] = sort(subMS,'descend');
for e = 1:100%size(dataStack_T,3)
    
    
    [dataOut] = getFramePackage(mC(:,ssidx(e)),mE,mU,szT(1:3));
    
    rawImage = [];
    simImage = [];
    for s = 1:size(dataStack_T,4)
        rawImage = cat(2,rawImage,bindVec(dataStack_T(:,:,ssidx(e),s)));
        simImage = cat(2,simImage,bindVec(mS(:,:,ssidx(e),s)));
    end
    imshow([simImage;rawImage],[]);
    
    
    hold on
    plot(mA1(ssidx(e),1),mA1(ssidx(e),2),'r.')
    plot(mA2(1,1,ssidx(e)),mA2(1,2,ssidx(e)),'g.')
    plot(mA2(2,1,ssidx(e)),mA2(2,2,ssidx(e)),'g.')
    hold off
    waitforbuttonpress
    drawnow
    
    
end
%% measure the raw data patches
close all
binTH = .1;
disp = 200;
NP = [];
MS = [];
CM = [];
EP_matrix = [];
domainArea = [];
frameImage = zeros(2*DMSZ+1,2*DMSZ+1);
for e = 1:4
    frameImage(:,1) = 1;
    frameImage = imrotate(frameImage,90); 
end

for e = 1:size(dataStack_T,3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skeleton measurements
    tmp = dataStack_T(:,:,e,4);
    tmp = logical(tmp > binTH);
    tmp = bwmorph(tmp,'skel',inf);
    tmp = bwmorph(tmp,'spur',3);
    ep = bwmorph(tmp,'endpoints');
    [e2,e1] = find(ep);
    eidx = find(ep);
    outOver1 = tmp;
    EP = [e1 e2];
    kpidx = frameImage(eidx) == 1;
    EP = EP(kpidx,:);
    
    
    NP(e) = size(EP,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % domain measurements
    tmp = logical(dataStack_T(:,:,e,2) > binTH);
    [d2,d1] = find(tmp);
    CM(e,:) = [mean(d1) mean(d2)];
    CM(e,:) = [21 21];
    domainArea(e) = sum(tmp(:));
    tmp1 = dataStack_T(:,:,e,end);
    tmp2 = dataStack_T(:,:,e,end-1);
    qX = [];
    qX(:,1) = ba_interp2(tmp1,EP(:,1),EP(:,2));
    qX(:,2) = ba_interp2(tmp2,EP(:,1),EP(:,2));
    qX_NOR = sum(qX.*qX,2).^-.5;
    qX = bsxfun(@times,qX,qX_NOR);
    
    qM = [];
    qM(:,1) = ba_interp2(tmp1,CM(e,1),CM(e,2));
    qM(:,2) = ba_interp2(tmp2,CM(e,1),CM(e,2));
    qM_NOR = sum(qM.*qM,2).^-.5;
    qM = bsxfun(@times,qM,qM_NOR);
    
    
    
    
    
    
    if NP(e) == 2
        vec1 = CM(e,:) - EP(1,:);
        vec2 = -(CM(e,:) - EP(2,:));
        vec1 = vec1 / norm(vec1);
        vec2 = vec2 / norm(vec2);
        MS(e) = vec1*vec2';
        EP_matrix = cat(3,EP_matrix,EP);
    elseif NP(e) == 1
        vec1 = CM(e,:) - EP(1,:);
        MS(e) = 1;
        EP_matrix = cat(3,EP_matrix,[EP;zeros(1,2)]);
    else 
        MS(e) = 0;
        EP_matrix = cat(3,EP_matrix,zeros(2));
    end
    
    
    if e < disp
        out = flattenMaskOverlay(dataStack_T(:,:,e,1),outOver1,1,'r');
        imshow(out,[]);
        hold on
        quiver(EP(:,1),EP(:,2),qX(:,1),qX(:,2),'m');
        quiver(EP(:,1),EP(:,2),-qX(:,1),-qX(:,2),'m');
        quiver(CM(e,1),CM(e,2),qM(1),qM(2),10,'m');
        quiver(CM(e,1),CM(e,2),-qM(1),-qM(2),10,'m');
          plot(EP(:,1),EP(:,2),'g.');
        plot(CM(e,1),CM(e,2),'b.');
        hold off
        drawnow
        axis([-5 50 -5 50])
        waitforbuttonpress
    end
   
    
   
    %waitforbuttonpress
end
