FilePath = '/mnt/snapper/nate/Madison 19 research/';
FileList = {};
FileExt = {'tif'};
verbose = 1;
FileList = sdig(FilePath,FileList,FileExt,verbose);
%% conformal map2
for e = 1:10
close all
D = generateSdomain(DMSZ,true);
pause(.5)
end
%%
close all
vxMax = 1;
ax = 2;
vXmin = 0;
breakPointTX = 0;
VX = @(y)(vxMax*(1+exp(-ax*(y-breakPointTX))).^-1) + vXmin;
vyMax = 1;
ay = .1;
vYmin = 0;
VY = @(y)(vyMax*(1+exp(-ay*y)).^-1) + vYmin;
VY = @(y)vyMax*ones(size(y));
initX = [0 0];

[x1,x2] = ndgrid(linspace(-DMSZ,DMSZ,2*DMSZ+1),linspace(-DMSZ,DMSZ,2*DMSZ+1));

nx2 = cumsum(VX(x1),1)+x2;
nx1 = x1(1,1)+cumsum(VY(x1),1);



for s = 1:size(nx1,1)
    plot(nx2(s,:),nx1(s,:),'r')
    hold all
end

for s = 1:size(nx1,1)
    plot(nx2(:,s),nx1(:,s),'b')
    hold all
end

%% conformal map1
close all
maxCompressValue = .1;
k = 1;
minCompressValue = .1;
alpha = @(y)(maxCompressValue*(1+exp(k*y)).^-1) + minV;

maxKurvature = .1;
kc = .1;
minKurvature = -.01;
c = @(y)(maxKurvature*(1+exp(-kc*y)).^-1) + minKurvature;
y = linspace(-10,10,100);
plot(y,alpha(y));
waitforbuttonpress
plot(y,c(y));
waitforbuttonpress
[x1,x2] = ndgrid(linspace(-DMSZ,DMSZ,2*DMSZ+1),linspace(-DMSZ,DMSZ,2*DMSZ+1));

close all

nx2 = alpha(x1).*x2;
nx1 = c(x1).*x2.^2 + x1;



for s = 1:size(nx1,1)
    plot(nx2(s,:),nx1(s,:),'r')
    hold all
end

for s = 1:size(nx1,1)
    plot(nx2(:,s),nx1(:,s),'b')
    hold all
end
%%
close all
dataStack = [];
dataStack_bin = [];
dataStack_dT = [];
samCNT = 1;
DMSZ = 20;
BUFFER = 20;
CH = 100;
perCH = 50;
%CH = 100;
[n1,n2] = ndgrid(-(DMSZ+BUFFER):(DMSZ+BUFFER),-(DMSZ+BUFFER):(DMSZ+BUFFER));
sDomain = [n2(:) n1(:) ones(size(n1(:)))];
%sDomain = [nx1(:) nx2(:) ones(size(n1(:)))];
dataStack = zeros([size(n1) CH]);
dataStack_bin = zeros([size(n1) CH]);
dataStack_dT = zeros([size(n1) CH]);
dataStack_skel = zeros([size(n1) CH]);
dataStack_nor1 = zeros([size(n1) CH]);
dataStack_nor2 = zeros([size(n1) CH]);
s = 1;
%for s = 1:numel(FileList)
stopFlag = false;
disp = 10;
while ~stopFlag
    try
        
        fprintf(['Loading images....START\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load the images
        out = [];
        I = [];
        for e = 1:numel(FileList{s})
            I(:,:,e) = imread(FileList{s}{e});
        end
        I = uint8(I);
        % load the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Loading images....END\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % thresholds for binary object removal
        % area threshold
        areaThreshold = 20;
        eccentricityThreshold = .93;
        motionAreaThreshold = 10;
        dilateValue = 5;
        spurValue = 5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Making masks....START\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the back ground image
        mskF = [];
        for t = 1:size(I,3)
            fprintf(['(' num2str(t) ')-->start\n']);tic
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the foreground
            tmp = double(I(:,:,t));
            tmp = imfilter(tmp,fspecial('average',61),'replicate');
            tmp = imdilate(tmp,strel('disk',71));
            tmp = imfilter(tmp,fspecial('average',61),'replicate');
            foreGround(:,:,t) = bindVec(double(I(:,:,t)) - tmp);
            msk1 = double(I(:,:,t)/255) < graythresh(double(I(:,:,t)/255));
            msk2 = double(foreGround(:,:,t)) < graythresh(foreGround(:,:,t));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % fill holes
            msk2 = imclose(msk2,strel('disk',dilateValue,0));
            msk2 = imfill(msk2,'holes');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clear sides
            bottom = msk2(end,:);
            msk2(end,:) = 0;
            msk2 = imclearborder(msk2);
            msk2(end,:) = bottom;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % apply area threshold
            msk2 = bwareaopen(msk2,areaThreshold);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clear sides
            out = flattenMaskOverlay(double(I(:,:,t))/255,msk1,.2,'r');
            out = flattenMaskOverlay(out,msk2,.2,'b');
            R = regionprops(msk2,'Perimeter','MajorAxisLength','MinorAxisLength',...
                                'Centroid','Area','Eccentricity','PixelIdxList');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % apply area threshold
            kidx1 = find([R.Eccentricity] > eccentricityThreshold);
            mskT = zeros(size(msk2));
            for obj = 1:numel(kidx1)
                mskT(R(kidx1(obj)).PixelIdxList) = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % store the mask
            mskF(:,:,t) = mskT;
            fprintf(['(' num2str(t) ')-->end @' num2str(toc) '\n']);
        end
        fprintf(['Making masks....END\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                R = regionprops(logical(objFilter),'BoundingBox');
                dT = double(bwdist(objFilter));
                dT = imfilter(dT,fspecial('gaussian',[21 21],4),'replicate');
                dK = imcrop(dT,R.BoundingBox);
                para.scales.value = 1;
                para.resize.value = 1;
                [K,lam] = surKur(dK,para);


                [dX,dY] = gradient(dT);

                dX = zeros(size(dX));
                dX(R.BoundingBox(2):(R.BoundingBox(2)+size(lam,1)-1),...
                    R.BoundingBox(1):(R.BoundingBox(1)+size(lam,2)-1)) = lam(:,:,1);
                dY = zeros(size(dY));
                dY(R.BoundingBox(2):(R.BoundingBox(2)+size(lam,1)-1),...
                    R.BoundingBox(1):(R.BoundingBox(1)+size(lam,2)-1)) = lam(:,:,2);
                


                
                
                
                
                
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
                        aff = [[nor';0] [tan';0] [m';1]];
                        iaff = inv(aff);










                        %[sDomain,dDomain] = generateSdomain(DMSZ,false);

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

                        % this works when the affine transformation is
                        % constant along the space
                        samD = [sam1(:) sam2(:) zeros(size(sam1(:)))]';
                        samD = (iaff*samD)';


                        %{
                        for trans = 1:size(dDomain,1)
                            tmpT = squeeze(dDomain(trans,:,:));
                            tmpT = [[tmpT zeros(size(tmpT,1),1)];[0 0 1]];
                            %tmpT = inv(tmpT);
                            %what
                            samD(:,trans) = tmpT*iaff*samD(:,trans);
                            %samD(:,trans) = iaff*tmpT*samD(:,trans);
                        end
                        samD = samD';
                        %}




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
            
           
            
            %drawnow
            %hold off
            %pause(.2)
            %waitforbuttonpress
        end
        
        
    catch ME
        getReport(ME)
        
    end
    
    s = s + 1;
end
%% here we go again
close all
DMSZ = 20;
dataStack_T = cat(4,dataStack,dataStack_bin,dataStack_dT,dataStack_skel,dataStack_nor2,dataStack_nor1);
mp = [(size(dataStack,1) - 1 )/2,(size(dataStack,2) - 1 )/2] + 1;



baseAffine = eye(3);
baseAffine(:,3) = [mp';1];


tr = 1;
    for rep = 1:100
    [D,transD,gD,sD] = generateSdomain(DMSZ,true);
    %D = sD;









    aD = (baseAffine*D')';
    imshow(dataStack_T(:,:,tr),[]);
    hold on
    plot(aD(:,1),aD(:,2),'.');
    gD = reshape(aD(:,1:2),size(gD));




    for s = 1:size(gD,1)
        plot(gD(s,:,1),gD(s,:,2),'r')
        hold all
    end
    for s = 1:size(gD,2)
        plot(gD(:,s,1),gD(:,s,2),'r')
        hold all
    end



    sam = squeeze(ba_interp2(dataStack_T(:,:,tr,:),aD(:,1),aD(:,2)));
    sam = reshape(sam,[2*DMSZ+1,2*DMSZ+1,size(dataStack_T,4)]);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skeleton measurements
    tmp = sam(:,:,4);
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
    NP = size(EP,1);




    tmp1 = sam(:,:,end);
    tmp2 = sam(:,:,end-1);


    % whole transformation
    for e = 1:size(transD,1)
        newT(e,:) = (squeeze(transD(e,:,:))*[tmp1(e);tmp2(e)])';
    end
    tmp1 = reshape(newT(:,1),size(tmp1));
    tmp2 = reshape(newT(:,2),size(tmp2));



    CM(e,:) = [21 21];



    qX = [];
    qX(:,1) = ba_interp2(tmp1,EP(:,1),EP(:,2));
    qX(:,2) = ba_interp2(tmp2,EP(:,1),EP(:,2));
    qX_NOR = sum(qX.*qX,2).^-.5;
    qX = bsxfun(@times,qX,qX_NOR);





    qM = [];
    qM(:,1) = ba_interp2(tmp1,CM(e,1),CM(e,2));
    qM(:,2) = ba_interp2(tmp2,CM(e,1),CM(e,2));
    %qM = -[1 .5];
    qM_NOR = sum(qM.*qM,2).^-.5;
    qM = bsxfun(@times,qM,qM_NOR);


    %qM = (inv(squeeze(transD(841,:,:)))*qM')';
    %qM = ((squeeze(transD(841,:,:)))*qM')';


    out = flattenMaskOverlay(sam(:,:,1),logical(tmp),1,'r');
    figure;
    imshow(out,[]);
    pause(.21)

    hold on
    quiver(EP(:,1),EP(:,2),qX(:,1),qX(:,2),'m');
    quiver(EP(:,1),EP(:,2),-qX(:,1),-qX(:,2),'m');


    quiver(CM(e,1),CM(e,2),qM(1),qM(2),10,'m');
    quiver(CM(e,1),CM(e,2),-qM(1),-qM(2),10,'m');

pause(1)
    close all
end

%% model all the data raw all
numC = 170;
dataStack_T = cat(4,dataStack,dataStack_bin,dataStack_dT,dataStack_skel,dataStack_nor2,dataStack_nor1);
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
binTH = .01;
disp = 40;
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
    

    % sim data
    tmp1S = mS(:,:,e,end);
    tmp2S = mS(:,:,e,end-1);
    qXS = [];
    qXS(:,1) = ba_interp2(tmp1S,EP(:,1),EP(:,2));
    qXS(:,2) = ba_interp2(tmp2S,EP(:,1),EP(:,2));
    qXS_NOR = sum(qXS.*qXS,2).^-.5;
    qXS = bsxfun(@times,qXS,qXS_NOR);
    
    
    
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

        % sim data
        quiver(EP(:,1),EP(:,2),qXS(:,1),qXS(:,2),10,'g');
        quiver(EP(:,1),EP(:,2),-qXS(:,1),-qXS(:,2),10,'g');

        plot(EP(:,1),EP(:,2),'g.');
        plot(CM(e,1),CM(e,2),'b.');
        hold off
        drawnow
        axis([-5 50 -5 50])
        waitforbuttonpress
    end
    %waitforbuttonpress
end
%% sub section the data
sidx = find(NP == 2);
%sidx = find(NP == 1);
numC = 20;
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
%% generate feature domains

sd = 11;
[w1,w2] = ndgrid(-sd:sd,-sd:sd);
sd = [w1(:) w2(:) ones(size(w1(:)))];
fdSZ = size(w1);
featureDomains{1} = sd;
featureDomains{2} = sd;
featureDomains{3} = sd;

%% scan data packages
al = [];
for e = 1:100
    [dataOut] = getDataPackage(mC(:,e),mE,mU,szT(1:3),featureDomains,fdSZ,true);
    al = [al dataOut.features(2)];
    %plot(al)
    %drawnow
end
%%
tr = 5;
mag = 2;
buildGrad(mag,mL.^.5,mC(:,tr),mE,mU,szT(1:3),featureDomains,fdSZ,true);
%% EXTRA PLAY
close all

subMS = MS(sidx);
[~,ssidx] = sort(subMS,'descend');

for e = 1:100%size(dataStack_T,3)
    
    mag = 2;
    buildGrad(mag,mL.^.5,mC(:,ssidx(e)),mE,mU,szT(1:3),featureDomains,fdSZ,true);
    
    
    [dataOut] = getDataPackage(mC(:,ssidx(e)),mE,mU,szT(1:3),featureDomains,fdSZ,true);
    
    
    [dataOut1] = getFramePackage(mC(:,ssidx(e)),mE,mU,szT(1:3),true);
    [dataOut2] = getFramePackage(mC(:,ssidx(e+1)),mE,mU,szT(1:3),true);
    
    
    
    %F = inv(dataOut2.inFrame)*dataOut2.cenFrame;
    %T21 = inv(dataOut1.cenFrame)*inv(dataOut2.inFrame)*dataOut2.cenFrame;
    %T21 = inv(dataOut2.inFrame)*dataOut2.cenFrame;
    %T = inv(inv(dataOut2.cenFrame)*dataOut2.inFrame);


    T = inv(dataOut2.inFrame)*dataOut2.cenFrame;
    F = dataOut1.cenFrame*T*inv(dataOut2.cenFrame)*dataOut2.inFrame;
    figure;hold on
    plotFrame(dataOut1.inFrame,10,{'r' 'r'});
    plotFrame(dataOut1.cenFrame,10,{'b' 'b'});
    plotFrame(dataOut1.outFrame,10,{'g' 'g'});
    plotFrame(F,10,{'m' 'm'});
    
     
     
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
