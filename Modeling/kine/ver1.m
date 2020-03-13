%% scan for data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FilePath = '/home/nate/RILPop/';
%FilePath = '/mnt/spaldingdata/Ashley/root_growth/RILpop/';
FileList = {};
FileExt = {'csv'};
%FileList = fdig(FilePath,FileList,FileExt,1);
FileList = gdig(FilePath,FileList,FileExt,1);
%% match up with the raw image set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FilePath = '/home/nate/RILPop/';
mFileList = {};
FileExt = {'mat'};
mFileList = fdig(FilePath,mFileList,FileExt,1);
%% raw image path
imgPath = '/mnt/spaldingdata/Ashley/root_growth/RILpop/';
%% try to auto click on the tip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the csv files based on this search
fidx1 = contains(FileList,'rawData');
fidx2 = contains(FileList,'RILpop');
fidx = fidx1 & fidx2;
FileList = FileList(fidx);
%% sample for training 1) tip finding 2) line tracing
% make grid for tip finding
samSZ = 100;
[s1,s2] = ndgrid(linspace(-samSZ,samSZ,2*samSZ+1),linspace(-samSZ,samSZ,2*samSZ+1));
samG = [s1(:) s2(:) ones(size(s1(:)))];
szG = size(s1);

% make grid for tracing
samSZ1 = 100;
[g1,g2] = ndgrid(linspace(-samSZ1,samSZ1,2*samSZ1+1),linspace(-samSZ1,samSZ1,2*samSZ1+1));
traceGrid = [g1(:) g2(:) ones(size(g1(:)))];
szG1 = size(g1);

% for each csv file

numToSample = 80;


tipSample = false;

if tipSample;X = [];Y = [];end
traceSample = false;
if traceSample;traceX = [];traceY = [];end

for e = 1:numToSample%numel(FileList)

    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the csv file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = csvread(FileList{e});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % string manipulation for the file path
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pth,nm,ext] = fileparts(FileList{e});
        didx = strfind(pth,'--');
        tnm = pth((didx(end-1)+2):(didx(end)-1));
        fidx = strfind(tnm,'_');
        l(e) = str2num(tnm((fidx(end)+1):end));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make mat file name
        % and load the mat data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        matName = [pth filesep 'rootTipObject.mat'];
        rootT = load(matName);
        P = rootT.d.position_field_midline;
        P = reshape(P.d,P.oS);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the image file name
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imgBase = strrep(pth,FilePath,'');
        jidx = strfind(imgBase,filesep);
        imgBase = imgBase(1:jidx(1)-1);
        imgBase = [imgPath imgBase filesep tnm];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dig for the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tFilePath = imgBase;
        tFileList = {};
        tFileExt = {'tif'};
        tFileList = fdig(tFilePath,tFileList,tFileExt,0);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % analysis of first frame
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mid = squeeze(P(:,1,:));
        mid = arcLength(mid,'arcLen');
        % read the image
        I = double(imread(tFileList{1}))/255;
        I = rgb2gray(I);
        %
        perT = .45;
        cI = imcomplement(I);
        rI = imresize(cI,.4);
        fI = imfilter(rI,fspecial('gaussian',[41 41],11),'replicate');
        fI2 = imopen(fI,strel('disk',41,0));
        %fI3 = imfilter(fI2,fspecial('gaussian',[41 41],11),'replicate');
        fI3 = imfilter(fI2,fspecial('disk',41),'replicate');
        BK = imresize(fI3,size(I));
        root = cI - BK;
        rootMask = root > graythresh(root)*perT;
        rootMask = bwlarge(rootMask);
        rootMask = imclose(rootMask,strel('disk',11,0));
        rootMask = imfill(rootMask,'holes');
        [midlineM,rootWidth(e),BOX,TIP{e},WHOLE{e},Domain] = extendMidline(255*double(rootMask.*I),mid');
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        midL = 1000;
        % make the mask - tmp
        tmpMask = WHOLE{e}~=0;
        % sum along x direction for width profile
        widthProfile = sum(tmpMask,2);
        % find where the profile starts
        fidx = find(widthProfile ~= 0);
        Emidline = squeeze(Domain(fidx(1):fidx(midL),(end-1)/2,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample the midline for tracing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if traceSample
            mSKIP = 20;
            fSKIP = 50;
            smoothV = 7;
            keepPer = 1/5;

            for f = 1:fSKIP:numel(tFileList)

                tan = [];nor = [];
                I = double(imread(tFileList{f}))/255;
                I = rgb2gray(I);
                mid = squeeze(P(:,f,:));

                mid = flip(Emidline,2);

                mid = imfilter(mid,ones(smoothV,1)/smoothV,'replicate');


                mid = arcLength(mid,'arcLen');

                tan = diff(mid,1,1);
                nor(:,1) = tan(:,2);
                nor(:,2) = -tan(:,1);
                pt = mid(1:(end-1),:);

                affine = vecStackToAffine(tan,nor,pt);
                amid = [mid ones(size(mid,1),1)];



                dL = 100;

                for p = 1:mSKIP:(size(affine,3)-dL)
                    tmpA = affine(:,:,p);
                    tmpMid = mtimesx(inv(tmpA),amid,'T')';
                    tmpD = mtimesx(tmpA,traceGrid,'T')';
                    tmpI = ba_interp2(I,tmpD(:,2),tmpD(:,1));
                    tmpI = reshape(tmpI,szG1);
                    hsz = hsize(tmpI);






                    tmpMid = tmpMid(p:(p+dL),:);
                    tmpMid(:,1:2) = bsxfun(@plus,tmpMid(:,1:2),hsz);


                    if rand(1) > keepPer
                        traceY = cat(3,traceY,tmpMid);
                        traceX = cat(3,traceX,tmpI);
                    end

                    if f == 1


                        imshow(I,[]);
                        hold on
                        plot(mid(:,2),mid(:,1),'b');
                        plot(mid(p:(p+dL),2),mid(p:(p+dL),1),'g')

                        imshow(tmpI,[]);
                        hold on
                        plot(tmpMid(:,2),tmpMid(:,1),'r');
                        hold off


                        drawnow
                    end




                end

                f

            end
        end



        if tipSample
            fr = 1;
            mid = squeeze(P(:,fr,:));
            mid = arcLength(mid,'arcLen');





            close all
            imshow(I,[]);hold on
            hold on;
            plot(mid(:,2),mid(:,1),'r');
            plot(mid(1,2),mid(1,1),'go');
            plot(Emidline(:,1),Emidline(:,2),'c');





            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample the tip only
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tp = [];ed = [];tan = [];nor = [];
            tipMask = zeros(size(I));
            tipMask(round(Emidline(1,2)),round(Emidline(1,1))) = 1;
            tipMask = imdilate(tipMask,strel('disk',5,0));
            [tp(:,1),tp(:,2)] = find(tipMask);


            E = edge(I,'Canny',[],15);
            E = bwlarge(E,5);
            E = imdilate(E,strel('disk',3,0));
            E = E & ~tipMask;

            [ed(:,1),ed(:,2)] = find(E);
            fI = imfilter(I,fspecial('gaussian',[31 31],7),'replicate');
            [d1,d2] = gradient(fI);
            dt = (d1.^2 + d2.^2).^.5;
            d1 = d1.*dt.^-1;
            d2 = d2.*dt.^-1;

            tan(:,2) = ba_interp2(d1,ed(:,2),ed(:,1));
            tan(:,1) = ba_interp2(d2,ed(:,2),ed(:,1));
            nor(:,1) = tan(:,2);
            nor(:,2) = -tan(:,1);

            %{
            imshow(I,[]);
            hold on
            quiver(ed(:,2),ed(:,1),tan(:,2),tan(:,1));
            quiver(ed(:,2),ed(:,1),nor(:,2),nor(:,1));
            %}
            affine = vecStackToAffine(tan,nor,ed);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample at near tip points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tipArea = bwdist(tipMask) < 50;
            %tipArea = tipArea & 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample at edge points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SAMSZ = 2000;
            SAMSZ = 300;
            samI = [];
            sampleSize = min(SAMSZ,size(affine,3));
            ridx = randperm(size(affine,3));
            affine = affine(:,:,ridx);
            affine = affine(:,:,1:sampleSize);
            samEDGE = zeros(size(samG,1),sampleSize);
            for p = 1:size(affine,3)
                tmpD = mtimesx(affine(:,:,p),samG,'T')';
                samEDGE(:,p) = ba_interp2(I,tmpD(:,2),tmpD(:,1));
            end
            samEDGE = reshape(samEDGE,[szG size(affine,3)]);
            vecEDGE = zeros(size(samEDGE,3),1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample at tip
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tan = [];nor = [];
            tan(:,2) = ba_interp2(d1,tp(:,2),tp(:,1));
            tan(:,1) = ba_interp2(d2,tp(:,2),tp(:,1));
            nor(:,1) = tan(:,2);
            nor(:,2) = -tan(:,1);
            affine = vecStackToAffine(tan,nor,tp);
            samTIP = zeros(size(samG,1),size(affine,3));
            for p = 1:size(affine,3)
                tmpD = mtimesx(affine(:,:,p),samG,'T')';
                samTIP(:,p) = ba_interp2(I,tmpD(:,2),tmpD(:,1));
            end
            samTIP = reshape(samTIP,[szG size(affine,3)]);

            vecTIP = ones(size(samTIP,3),1);

            X = cat(3,X,samEDGE,samTIP);
            Y = [Y;vecEDGE;vecTIP];
        end

        
        fprintf(['Done sampling:' num2str(e) ':' num2str(numToSample) '\n']);
    catch ME
        ME
    end
end
xSZ = size(X);
X = reshape(X,[xSZ(1:2) 1 xSZ(3)]);
%% try patten net and PCA
szX = size(X);
X = reshape(X,[prod(szX(1:3)) szX(end)]);
[xU,xE,xL] = PCA_FIT_FULL_Tws(X,5);
pX = PCA_REPROJ_T(X,xE,xU);
%% PATTERN NET - PCA
tipNet = patternnet([5]);
pY = full(ind2vec((Y+1)',2));
tipNet = train(tipNet,pX,pY);
%% construct the points to classify NOT WORKING AS WELL
xSZ = size(X);

options = trainingOptions('sgdm',...
    'InitialLearnRate',0.01,...
    'MaxEpochs',100,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel');

layers = [ ...
    imageInputLayer([xSZ(1:2) 1],'Normalization','none')
    convolution2dLayer(11,11)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)

    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

trainedNet = trainNetwork(X,categorical(Y),layers,options);
%% run the root tip finding method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 210%%numel(FileList)

    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the csv file
        d = csvread(FileList{e});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % string manipulation for the file path
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pth,nm,ext] = fileparts(FileList{e});
        didx = strfind(pth,'--');
        tnm = pth((didx(end-1)+2):(didx(end)-1));
        fidx = strfind(tnm,'_');
        l(e) = str2num(tnm((fidx(end)+1):end));

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make mat file name
        % and load the mat data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        matName = [pth filesep 'rootTipObject.mat'];
        rootT = load(matName);
        P = rootT.d.position_field_midline;
        P = reshape(P.d,P.oS);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the image file name
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imgBase = strrep(pth,FilePath,'');
        jidx = strfind(imgBase,filesep);
        imgBase = imgBase(1:jidx(1)-1);
        imgBase = [imgPath imgBase filesep tnm];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dig for the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tFilePath = imgBase;
        tFileList = {};
        tFileExt = {'tif'};
        tFileList = fdig(tFilePath,tFileList,tFileExt,0);
        
        
        tan = [];nor = [];ed = [];
        I = double(imread(tFileList{1}))/255;
        I = rgb2gray(I);

        E = edge(I,'Canny',[],15);
        E = bwlarge(E,3);
        E = imdilate(E,strel('disk',3,0));
        E = E & ~tipMask;
       
        [ed(:,1),ed(:,2)] = find(E);
        fI = imfilter(I,fspecial('gaussian',[31 31],7),'replicate');
        [d1,d2] = gradient(fI);
        dt = (d1.^2 + d2.^2).^.5;
        d1 = d1.*dt.^-1;
        d2 = d2.*dt.^-1;

        tan(:,2) = ba_interp2(d1,ed(:,2),ed(:,1));
        tan(:,1) = ba_interp2(d2,ed(:,2),ed(:,1));
        nor(:,1) = tan(:,2);
        nor(:,2) = -tan(:,1);

        %{
        imshow(I,[]);
        hold on
        quiver(ed(:,2),ed(:,1),tan(:,2),tan(:,1));
        quiver(ed(:,2),ed(:,1),nor(:,2),nor(:,1));
        %}
        affine = vecStackToAffine(tan,nor,ed);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample at edge points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        samI = [];
        samEDGE = zeros(size(samG,1),size(affine,3));
        for p = 1:size(affine,3)
            tmpD = mtimesx(affine(:,:,p),samG,'T')';
            samEDGE(:,p) = ba_interp2(I,tmpD(:,2),tmpD(:,1));
            p
            size(affine,3)
        end

        samEDGE = reshape(samEDGE,[szG size(affine,3)]);

        
        xSZ = size(samEDGE);
        samEDGE = reshape(samEDGE,[xSZ(1:2) 1 xSZ(3)]);


        samSZ = size(samEDGE);
        samEDGE = reshape(samEDGE,[prod(samSZ(1:3)) samSZ(end)]);
        pE = PCA_REPROJ_T(samEDGE,xE,xU);
        Ypre = tipNet(pE);
        Ypre = Ypre(2,:)';
        %Ypre = trainedNet.classify(samEDGE);

        rMask = zeros(size(I));
        ridx = sub2ind(size(rMask),ed(:,1),ed(:,2));
        rMask(ridx) = Ypre;

        sMask = imfilter(rMask,fspecial('gaussian',[41 41],5),'replicate');
        sMask = rMask;
        [m1,m2] = find(sMask == max(sMask(:)));


        close all
        imshow(I,[]);
        hold on;
        plot(ed(:,2),ed(:,1),'b.')
        plot(m2,m1,'ro')
        drawnow

    catch ME
        ME
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data - SIMPLE
% fit + root width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fidx1 = contains(FileList,'rawData');
%fidx2 = contains(FileList,'RILpop');
%fidx = fidx1 & fidx2;
%FileList = FileList(fidx);
func = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);

clear P
basePer = .1;
disp = false;
LEN = zeros(numel(FileList),1);
G = zeros(numel(FileList),4);
err = LEN;
clear TIP
close all
widthProfileS = [];
loadList = {};
failList = {};
rootWidthLoad = false;
LL = 500;
% for each csv file
for e = 1:numel(FileList)
    
    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the csv file
        d = csvread(FileList{e});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % string manipulation for the file path
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [pth,nm,ext] = fileparts(FileList{e});
        didx = strfind(pth,'--');
        tnm = pth((didx(end-1)+2):(didx(end)-1));
        fidx = strfind(tnm,'_');
        l(e) = str2num(tnm((fidx(end)+1):end));

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make mat file name
        % and load the mat data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        matName = [pth filesep 'rootTipObject.mat'];
        rootT = load(matName);
        P = rootT.d.position_field_midline;
        P = reshape(P.d,P.oS);
        
        if rootWidthLoad
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make the image file name
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imgBase = strrep(pth,FilePath,'');
            jidx = strfind(imgBase,filesep);
            imgBase = imgBase(1:jidx(1)-1);
            imgBase = [imgPath imgBase filesep tnm];


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % dig for the images
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tFilePath = imgBase;
            tFileList = {};
            tFileExt = {'tif'};
            tFileList = fdig(tFilePath,tFileList,tFileExt,0);
        end
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack whole movie
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DELTA = 1000;
        for f = 1:numel(tFileList)
            % read fth image
            I = double(imread(tFileList{f}))/255;
            mid = squeeze(P(:,f,:));
            mid = arcLength(mid,'arcLen');
            
            
            I = rgb2gray(I);
            perT = .45;
            cI = imcomplement(I);
            rI = imresize(cI,.4);
            fI = imfilter(rI,fspecial('gaussian',[41 41],11),'replicate');
            fI2 = imopen(fI,strel('disk',41,0));
            %fI3 = imfilter(fI2,fspecial('gaussian',[41 41],11),'replicate');
            fI3 = imfilter(fI2,fspecial('disk',41),'replicate');
            BK = imresize(fI3,size(I));
            root = cI - BK;
            rootMask = root > graythresh(root)*perT;
            rootMask = bwlarge(rootMask);
            rootMask = imclose(rootMask,strel('disk',11,0));
            rootMask = imfill(rootMask,'holes');

            rootMask = double(rootMask);
            maskStack{f} = rootMask;
            
            [~,~,~,~,tmp] = extendMidline(255*double(rootMask.*I),mid');

            
            tmpW = sum(tmp~=0,2);
            fidx = find(tmpW~=0);
            tmp = tmp(fidx(1):(fidx(1)+DELTA),:);
            imshow(tmp,[]);
            drawnow
            tmpS(:,:,f) = tmp;
        end
        
        close all
        for t = 1:size(tmpS,3)
            imshow(tmpS(:,:,t),[]);
            drawnow
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack whole movie
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        if rootWidthLoad
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % analysis of first frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mid = squeeze(P(:,1,:));
            mid = arcLength(mid,'arcLen');
            % read the image
            I = double(imread(tFileList{1}))/255;
            I = rgb2gray(I);
            %
            perT = .45;
            cI = imcomplement(I);
            rI = imresize(cI,.4);
            fI = imfilter(rI,fspecial('gaussian',[41 41],11),'replicate');
            fI2 = imopen(fI,strel('disk',41,0));
            %fI3 = imfilter(fI2,fspecial('gaussian',[41 41],11),'replicate');
            fI3 = imfilter(fI2,fspecial('disk',41),'replicate');
            BK = imresize(fI3,size(I));
            root = cI - BK;
            rootMask = root > graythresh(root)*perT;
            rootMask = bwlarge(rootMask);
            rootMask = imclose(rootMask,strel('disk',11,0));
            rootMask = imfill(rootMask,'holes');

            [midlineM,rootWidth(e),BOX,TIP{e},WHOLE{e}] = extendMidline(255*double(rootMask.*I),mid');


            % make the mask - tmp
            tmpMask = WHOLE{e}~=0;
            % sum along x direction for width profile
            widthProfile = sum(tmpMask,2);
            % find where the profile starts
            fidx = find(widthProfile ~= 0);
            % clip and store the profile
            widthProfileS(:,e) = widthProfile(fidx(1):fidx(LL));


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % vis for the root width mask
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%
            % vis for first frame
            widx = .5*mean(widthProfileS,2);
            wMax = 300;
            wi = zeros(LL,wMax);
            cp = (wMax-1)/2 + 1;
            rightH = [(1:size(wi,1))' cp + widx];
            leftH = [(1:size(wi,1))' cp - widx];
            rightH = flip(rightH,1);
            fullH = [leftH;rightH];
            wi = poly2mask(fullH(:,2),fullH(:,1),size(wi,1),size(wi,2));
            imshow(wi,[]);
            drawnow
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % analysis of first frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init guess for vel fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        initP = zeros(1,4);
        
        initP(4) = 1.2;
        initP(3) = .008;
        
        sd = sort(d(:,2),'descend');
        initP(2) = mean(sd(1:round(numel(sd)*basePer)));
        initP(1) = mean(d(:,1));
        toFit = @(P)func(d(:,1),P);
        delta = @(P)-log(mean(normpdf((toFit(P) - d(:,2)))));


        para = fminsearch(delta,initP);

        xpos = linspace(0,1500,1500);
        vel = func(xpos,para);

        
        G(e,:) = para;
        LEN(e) = l(e);
        err(e) = delta(para);
        
        if disp
            plot(d(:,1),d(:,2),'.','MarkerSize',1);
            hold on
            plot(xpos,vel,'r');
            axis([0 1500 0 5]);
            drawnow
            hold off
        end
        
        
        loadList{end+1} = FileList{e};

    catch ME
        widthProfileS(:,e) = NaN*ones(size(widthProfileS,1),1);
        G(e,:) = NaN*ones(1,4);
        LEN(e) = NaN;
        err(e) = NaN;
        failList{end+1} = FileList{e};
    end
    e
end
%% remove nan - no need to remove from loadList
rmidx = any(isnan(G),2);
G(rmidx,:) = [];
LEN(rmidx) = [];
err(rmidx) = [];
%%
widthProfileS(:,rmidx) = [];
%% outlier remove
TF = [];
for e = 1:size(G,2)
    TF(:,e) = isoutlier(G(:,e));
end
rmidx = any(TF,2);
G(rmidx,:) = [];
LEN(rmidx) = [];
err(rmidx) = [];
%%
widthProfileS(:,rmidx) = [];
%% remove outlier based on width profile
[~,widC,widU,widE,widL,widERR,widLAM] = PCA_FIT_FULL(widthProfileS',2);
TF = [];
for e = 1:size(widC,2)
    TF(:,e) = isoutlier(widC(:,e));
end
rmidx = any(TF,2);
G(rmidx,:) = [];
LEN(rmidx) = [];
err(rmidx) = [];
widthProfileS(:,rmidx) = [];
%% cell length
X = linspace(0,2500,2500);
regr = regrFunc(X,G(1,:));
plot(X,regr)
[i1,i2] = ndgrid(linspace(.2*10^-3,1*10^-3,25),linspace(.1,1,25));
initP = genCellLength(velFunc,regrFunc,X,G(tr,:),initP,dt,TAU,2.^-cellD);
%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert
conFPS = 30^-1;
conSPH = 60*60;
conPPM = 1463^-1;
unitsG = G;
unitsG(:,2) = G(:,2)*conFPS*conSPH*conPPM;
unitsG(:,1) = G(:,1)*conPPM;
units_to_org = @(p)bsxfun(@times,p,[conPPM^-1 (conFPS*conSPH*conPPM)^-1 1 1]);
%% sweep original parameters
swG = std(G,1,1);
uG = mean(G,1);
N = 5;
mag = 1.5;
X = linspace(0,2500,2500);
velFunc = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);
regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;
LABEL = {'x0','vf','k','n','len'};
close all
for p = 1:size(uG,2)
    L = linspace(-mag*swG(p),mag*swG(p),N);
    for step = 1:N
        tmpP = uG;
        tmpP(p) = tmpP(p) + L(step);
        Yvel = velFunc(X,tmpP);
        Yregr = regrFunc(X,tmpP);
        yyaxis('left');
        plot(X,Yvel,'k-');
        xlabel('Position (px)')
        ylabel('Velocity (px/fr)')
        hold on
        yyaxis('right');
        plot(X,Yregr,'r-');
        ylabel('REGR (%/fr)')
        hold on
    end
    title(LABEL{p})
    hold off
    waitforbuttonpress
    close all
end
%% pairwise plots - for unit converted data
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(unitsG,2,LABEL);
%% make PCA traces on data
[zunitsG,zuG,zsG] = zscore(unitsG);
z_to_org = @(p)bsxfun(@plus,bsxfun(@times,p,zsG),zuG);

[unit_latS,unit_latC,unit_latU,unit_latE,unit_latL,unit_latERR,lunit_atLAM] = PCA_FIT_FULL(zunitsG,3);
%unit_latE = bsxfun(@times,unit_latE',zsG)';
%unit_latU = mean(unitsG);
mag = 5;
unitCurves = {};
for e = 1:size(unit_latE,2)
    lunit_atLAM(e,e)
    
    
    dL = linspace(-mag*lunit_atLAM(e,e).^-.5,mag*lunit_atLAM(e,e).^-.5,7);
    dL = bsxfun(@plus,(unit_latE(:,e)*dL)',unit_latU);
    dL = z_to_org(dL);
    unitCurves{e} = dL;
    
    
end
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(unitsG,2,LABEL,{unitCurves});
%% make new curves for tracing
close all
[Z,tU,tS] = zscore(G);
[~,tC,tU,tE,tL,tERR,tLAM] = PCA_FIT_FULL(Z,1);
tE = (bsxfun(@times,tE',tS)');
tE = tE / norm(tE);
%%
nG = bsxfun(@minus,G,mean(G,1));
s = std(nG*tE,1,1);
tpc = linspace(-s(1),s(1),2);
tpc = bsxfun(@plus,(tE(:,1)*tpc)',mean(G,1));

P = nchoosek((1:size(G,2)),2);
for e = 1:size(P,1)
    tmpX = G(:,P(e,:));
    tmpV = tE(P(e,:)',1);
    tmpG = G(:,P(e,:));
    tmpV = tmpV / norm(tmpV);
    
    ntmpG = bsxfun(@minus,tmpG,mean(tmpG,1));
    
    
    s = std(ntmpG*tmpV);
    tpc = linspace(-s,s,3);
    tpc = bsxfun(@plus,(tmpV*tpc)',mean(tmpG,1));
    
    plot(tmpG(:,1),tmpG(:,2),'k.');
    hold on
    plot(tpc(:,1),tpc(:,2),'r');
    waitforbuttonpress
    close all
end

%%

for c = 1:numel(unitCurves)
    

    
    
    
    for p = 1:size(unitCurves{c},1)
        
        tmpP = unitCurves{c}(p,:);
        
        vel = velFunc(X,tmpP);
        regr = regrFunc(X,tmpP);
        
        vel = vel*conFPS*conSPH*conPPM;
        regr = 100*regr*conFPS*conSPH;
        pX = X*conPPM;
        
        yyaxis('left');
        hold on
        plot(pX,vel,'k-');
       
        xlabel('Position (mm)');
        ylabel('Velocity (mm/hr)');
        
        yyaxis('right');
        hold on
        plot(pX,regr,[CL{c} '-'])
        ylabel('REGR (%/hr)');
    end
    
    
    P = nchoosek((1:size(G,2)),2);
   
   
    for e = 1:size(P,1)
        figure;
        tX = G(:,P(e,:));
        plot(tX(:,1),tX(:,2),'k.');
        xlabel(LABEL{P(e,1)});
        ylabel(LABEL{P(e,2)});
        hold on
        
        for c = 1:numel(unitCurves)
            cur = unitCurves{c};
            tmpC = cur(:,[P(e,1),P(e,2)]);
            plot(tmpC(:,1),tmpC(:,2),CL{c});
        end
        
        xlabel(LABEL{P(e,1)});
        ylabel(LABEL{P(e,2)});
        waitforbuttonpress
    end
    
    
    hold off
    waitforbuttonpress
    close all

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%% pca to reduce from 4-->3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XalongRoot = linspace(0,3000,100);
[zG,zU,zS] = zscore(G);
[latS,latC,latU,latE,latL,latERR,latLAM] = PCA_FIT_FULL(G,3);
aff = zeros(1,size(latE,2)+1);
aff(end) = 1;
affineM = [[latE,latU'];aff];
affineMz = [[latE*latLAM.^.5,latU'];aff];
velFunc = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);
regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;
bk_lat_to_org = @(p)(mtimesx(affineM,p'))';
bk_lat_to_org_simple = @(p)(mtimesx(latE,p')'+latU);
re_org_to_lat_simple = @(p)(mtimesx(latE,'T',(p-latU)'))';
re_org_to_lat = @(p)(mtimesx(pinv(affineM),p'))';
bk_nlat_to_org = @(p)(mtimesx(affineMz,p'))';
bk_org_to_nlat = @(x)(mtimesx(inv(affineMz),p'))';
z_to_org = @(p)bsxfun(@plus,bsxfun(@times,p,zS),zU);
org_to_z = @(p)bsxfun(@times,bsxfun(@minus,p,zU),zS.^-1);
orgVec_to_z = @(p)bsxfun(@times,p,zS.^-1);
zVec_to_org = @(p)bsxfun(@times,p,zS);
%% map the data points to the curvilinear space from the parameter space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear p2m
tic;[mG] = p2m(G,XalongRoot);toc
tic;[mGn] = p2mn(G);toc
% NOTE: which conversion to use !!!!!
[zG,curU,curSIG] = zscore(mG);
% eigen vectors in the curvilinear space
[curS,curC,curU,curE,curL,curERR,curLAM] = PCA_FIT_FULL(mG,4);
% look at the z-score eigenvects
[zcurS,zcurC,zcurU,zcurE,zcurL,zcurERR,zcurLAM] = PCA_FIT_FULL(zG,4);
zcurE = bsxfun(@times,zcurE',curSIG)';
%% width profile analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
[~,widC,widU,widE,widL,widERR,widLAM] = PCA_FIT_FULL(widthProfileS',2);
plot(curC(:,1),curC(:,2),'.');
waitforbuttonpress
plot(widC(:,1),curC(:,1),'.');
waitforbuttonpress
figure;
plot(widC(:,1),mGn(:,1),'.');
figure;
plot(widC(:,1),G(:,1),'.');
plot(widC(:,1),LEN,'.');
%% CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = mGn;
Y = widC;
[X,xU,xSIG] = zscore(X);
[Y,yU,ySIG] = zscore(Y);
[A,B,r,U,V] = canoncorr(X,Y);
A = diag(xSIG)*A;
B = diag(ySIG)*B;
%% direct sum of kinematics + width profile - PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
X = [mGn,widC];
[~,dpC,dpU,dpE,dpL,dpERR,dpLAM] = PCA_FIT_FULL(X,3);
a = widC*dpE(5:6,1);
b = mGn*dpE(1:4,1);
plot(a,b,'.');
%% direct sum of kinematics + width profile - PCA - zscore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
X = [mGn,widC];
[X,xU,xSIG] = zscore(X);
[~,dpC,dpU,dpE,dpL,dpERR,dpLAM] = PCA_FIT_FULL(X,3);
dpE = diag(xSIG)*dpE;
a = widC*dpE(5:6,1);
b = mGn*dpE(1:4,1);
plot(a,b,'.');
%% width profile prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = G;
%X = mGn;
%X = [G mGn];
%X = curC;
%X = mGn;
Y = widC(:,2);
rq1 = [ones(size(X,1),1) X]\widC;
b = robustfit(X,Y);
widPre = [ones(size(X,1),1) X]*rq1;
widPre2 = [ones(size(X,1),1) X]*b;
close all
%plot(Y,widPre(:,1),'.');
hold on
plot(Y,widPre2(:,1),'r.');
plot([-800 800],[-800 800],'g')
v1 = b(2:end);
v1 = v1 / norm(v1);
v2 = curE(:,4);
v2 = v2 / norm(v2);
v1
v2
%% back predict peak (pc4) from width
%%%%%%%%%%%%%%%%%%%%
close all
X = widC(:,1:2);
Y = curC(:,4);
Y = mGn(:,3);
Y = zscore(Y);
X = zscore(X);
b = robustfit(X,Y);
%b = X\Y;
preY = b(1) + X*b(2:end);
%preY = X*b;
plot(Y,preY,'.')
%% width profile vs length
%%%%%%%%%%%%%%%%%%%%
close all
plot(LEN,widC(:,1),'.')
%% width profile sweep
%%%%%%%%%%%%%%%%%%%%
close all
[sweepD] = sweepPCA(widC,widE,widU,diag(widLAM).^.5,[1:2],5);
for c = 1:size(sweepD,1)
    for step = 1:size(sweepD,2)
        widthC = squeeze(sweepD(c,step,:));

        widx = .5*widthC;
        wMax = 300;
        wi = zeros(1000,wMax);
        cp = (wMax-1)/2 + 1;
        rightH = [(1:size(wi,1))' cp + widx];
        leftH = [(1:size(wi,1))' cp - widx];
        rightH = flip(rightH,1);
        fullH = [leftH;rightH];
        wi = poly2mask(fullH(:,2),fullH(:,1),size(wi,1),size(wi,2));
        imshow(wi,[]);
        drawnow
        waitforbuttonpress
    end
end
%% look at measurements slow vs fast
for e = 1:size(mGn,2)
    plot(mGn(:,e),mG(:,e),'.')
    waitforbuttonpress
end
close all
%% assign mGn to mG
%mG = mGn;
%% pairwise plots - save data to location for power point
% FIG #1
toSave = '/home/nate/forA/';
preFix = 'parameterSpace_pairwise';
toWait = false;
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(G,2,LABEL,[],figure,[],toSave,preFix,toWait);
close all
CON = {'peakPosition','peak','Width','asym'};
preFix = 'measureSpace_pairwise';
pairWisePlots(mGn,2,CON,[],figure,[],toSave,preFix,toWait);
close all
%% try to inverse map through the sub space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear P
velFunc_p = @(X,P)velFunc(X,bk_lat_to_org(P));
regrFunc_p = @(X,P)regrFunc(X,bk_lat_to_org(P));
regrPer = .8;
% this line was converted so that number of points along the domain do not
% influence the results
%toMetricSpace_p = @(P)extractREGRmetrics(@(P)regrFunc_p(XalongRoot,P),P,regrPer);
toMetricSpace_p = @(P)p2m(bk_lat_to_org(P),XalongRoot);
toMetricSpace_p2 = @(P)p2mn(bk_lat_to_org(P));
toMetricSpace_pn = @(P)p2mn(bk_nlat_to_org(P));
toMetricSpace_p3 = @(P)p2mn(bk_lat_to_org_simple(P));
toMetricSpace_raw = @(P)p2mn(P);
toMetricSpace_ZT = @(P)p2mn(z_to_org(P));
%% find inf the location of the mean
initP = [0;0;0;1];
norCur = @(x)bsxfun(@times,bsxfun(@minus,x,curU),curSIG.^-1); % one is added for affine mapping
uean_AtMetric = norCur(curU);
[uean_AtPara] = m2p(uean_AtMetric,@(v)norm(v),initP,XalongRoot,bk_nlat_to_org,norCur);
uean_AtPara = bk_nlat_to_org(uean_AtPara);
uean_AtPara(end) = [];
%{
%% trace vector curves - from CCA, direct-sum PCA, or PCA on zscores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: if the resolution of the p2mn is too low the Jacobian does not work
traceC = {};
traceVecSet = zcurE;
traceVecSet = curE;
%traceVecSet = A(:,1);
%tarceVecSet = dpE(1:4,1:3);
nTrace = 30;
dv = traceVecSet(:,1);
alpha = .1;
traceInLat1_pos = org_to_z(uean_AtPara);
[traceInLat1_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_pos,dv,nTrace,alpha);
traceInLat1_neg = org_to_z(uean_AtPara);
[traceInLat1_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_neg,-dv,nTrace,alpha);
traceC{1} = z_to_org([flip(traceInLat1_pos,1);traceInLat1_neg]);
%
dv = traceVecSet(:,2);
alpha = .1;
traceInLat2_pos = org_to_z(uean_AtPara);
[traceInLat2_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat2_pos,dv,nTrace,alpha);
traceInLat2_neg = org_to_z(uean_AtPara);
[traceInLat2_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat2_neg,-dv,nTrace,alpha);
traceC{2} = z_to_org([flip(traceInLat2_pos,1);traceInLat2_neg]);
%
dv = traceVecSet(:,3);
alpha = .1;
traceInLat3_pos = org_to_z(uean_AtPara);
[traceInLat3_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat3_pos,dv,nTrace,alpha);
traceInLat3_neg = org_to_z(uean_AtPara);
[traceInLat3_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat3_neg,-dv,nTrace,alpha);
traceC{3} = z_to_org([flip(traceInLat3_pos,1);traceInLat3_neg]);
%}
%% trace vector curves - from CCA, direct-sum PCA, or PCA on zscores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: if the resolution of the p2mn is too low the Jacobian does not work
traceC = {};
traceVecSet = zcurE;
%traceVecSet = curE;
%traceVecSet = A(:,1);
%tarceVecSet = dpE(1:4,1:3);
nTrace = 30;
for e = 1:size(traceVecSet,2)
    dv = traceVecSet(:,e);
    alpha = .1;
    traceInLat1_pos = org_to_z(uean_AtPara);
    [traceInLat1_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_pos,dv,nTrace,alpha);
    traceInLat1_neg = org_to_z(uean_AtPara);
    [traceInLat1_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_neg,-dv,nTrace,alpha);
    traceC{e} = z_to_org([flip(traceInLat1_pos,1);traceInLat1_neg]);
end
ztraceC = traceC;
%% trace vector curves - from CCA, direct-sum PCA, or PCA on zscores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: if the resolution of the p2mn is too low the Jacobian does not work
traceC = {};
%traceVecSet = zcurE;
traceVecSet = curE;
%traceVecSet = A(:,1);
%tarceVecSet = dpE(1:4,1:3);
nTrace = 30;
for e = 1:size(traceVecSet,2)
    dv = traceVecSet(:,e);
    alpha = .1;
    traceInLat1_pos = org_to_z(uean_AtPara);
    [traceInLat1_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_pos,dv,nTrace,alpha);
    traceInLat1_neg = org_to_z(uean_AtPara);
    [traceInLat1_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_neg,-dv,nTrace,alpha);
    traceC{e} = z_to_org([flip(traceInLat1_pos,1);traceInLat1_neg]);
end
ntraceC = traceC;
%% start gravitropism here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% place the left and right endpoints
% try to first place the endpoints along the eigencurves
% later move to a different perspective


X = linspace(0,2500,2500);    
alpha = .01;
mag = linspace(1,2,10);



traceVecSet = zcurE(:,1);

nTrace = 100;
traceGravi = {};
for e = 1:size(traceVecSet,2)
    dv = traceVecSet(:,e);

    traceInLat1_pos = org_to_z(uean_AtPara);
    [traceInLat1_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_pos,dv,nTrace,alpha);
    traceInLat1_neg = org_to_z(uean_AtPara);
    [traceInLat1_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_neg,-dv,nTrace,alpha);
    traceInLat1_neg(1,:) = [];
    traceGravi{e} = z_to_org([flip(traceInLat1_pos,1);traceInLat1_neg]);
end
traceGravi = traceGravi{1};
%%
close all
fL = 7;
curveT = imfilter(traceGravi,ones(fL,1)/fL,'replicate');
dL = gradient(curveT')';
dL = sum(dL.^2,2).^.5;
dL = cumsum(dL);
dL = dL - dL((end-1)/2);
plot(dL)

%%
close all
maxL = 10.5;
PP = linspace(0,pi,1000);
iC = interp1(dL,curveT,0);

midlineP = [X',zeros(size(X')),zeros(size(X'))];
tipAngle = 0;
dW = 90;
for tm = 1%:numel(PP)
    aL(tm) = sin(PP(tm));
    
    iL = maxL*sin(PP(tm));
    iR = -iL;
    iL = interp1(dL,curveT,iL);
    iR = interp1(dL,curveT,iR);
    
    
    tmp = [iL;iC;iR];
    
    tmp(:,2) = tmp(2,2);
    
    newV = domainFunctionAtPVect(X,tmp,velFunc);
    %plot(newV)
    %waitforbuttonpress
    vort = .5*gradient(newV,dW);
    
    dA = vort(:,2)*max(X)^-1;
    dA = cumsum(dA);
    dA = dA / 10000;
    for q = 1:10%(size(midlineP,1)-1)
        p = q + 1;
        dp = midlineP(p,3) + dA(p);
        
        tmpM = midlineP(p:end,1:2);
        cp = tmpM(1,1:2);
        tmpM = bsxfun(@minus,tmpM,cp);
        
        tmpR = sum(tmpM.*tmpM,2).^.5;
        
        dp = dp*ones(size(tmpR));
        dx = tmpR.*cos(dp);
        dy = tmpR.*sin(dp);
        
        %delta = [dx dy] - tmpM;
        %bsxfun(@plus,tmpM,cp);
        %midlineP(p,:) = midlineP(p,:) + [dx dy dp];
        
        %midlineP(p:end,:) = midlineP(p:end,:) + [delta dp];
        midlineP(p:end,3) = midlineP(p:end,3) + [dp];
        midlineP(p:end,1:2) = [dx dy];
        if any(isnan(midlineP(:)))
            break
        end
    end
    %waitforbuttonpress
    plot(midlineP(:,1),midlineP(:,2),'.')
    drawnow
    %waitforbuttonpress
    
    dA = max(X)^-1*(cumsum(vort(:,2))*180/pi);
    
    
    tipAngle(tm+1) = tipAngle(tm) + dA(end);
    
    %plot(tipAngle)
    %drawnow
    %%
    %{
    initP = traceGravi;
    close all
    options = optimset('Display','Iter','TolFun',10^-1000,'TolX',10^-1000,'MaxFunEvals',5000,'MaxIter',2000); 

    funcky = @(dP)graviContrast(X,initP,dP,velFunc);
    dP1 = fminsearch(funcky,zeros(1,size(initP,2)),options);
    [d,newP] = funcky(dP1);
    %traceGravi(2,:) = traceGravi(2,:) + dP;
    newV = domainFunctionAtPVect(X,newP,velFunc);
    vort = .5*gradient(newV);
    %}
   
    %initP = newP;
end

%%
%% end gravitropism here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make traces along the u-space vecs
nTrace = 30;
DV = eye(size(mG,2));
for v = 1:size(DV,2)
    dv = DV(:,v);
    alpha = .1;
    tmpTrace = org_to_z(uean_AtPara);
    [tmpTracePos] = traceVec_ver2(toMetricSpace_ZT,tmpTrace,dv,nTrace,alpha);
    [tmpTraceNeg] = traceVec_ver2(toMetricSpace_ZT,tmpTrace,-dv,nTrace,alpha);
    traceI{v} = z_to_org([flip(tmpTracePos,1);tmpTraceNeg]);
end
%% make traces of eigenVectors
% FIG #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make curves from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st = diag(curLAM).^.5;
vec = curE;
traceC = [];
for e = 1:size(vec,2)
    l = 30*linspace(-st(e),st(e),2);
    v = vec(:,e);
    v = (v*l)';
    v = bsxfun(@plus,v,curU);
    traceC{e} = v;
end
% save figures
preFix = 'eigenVectors_in_measureSpace';
CON = {'peakPosition','peakValue','Width','asym'};
pairWisePlots(mGn,2,CON,{traceC},figure,[],toSave,preFix,toWait);
close all
%% make traces of eigenVectors - zscores
% FIG #3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make eigen curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
st = diag(zcurLAM).^.5;
vec = zcurE;
traceCZ = [];
for e = 1:size(vec,2)
    l = 30*linspace(-st(e),st(e),2);
    v = vec(:,e);
    v = (v*l)';
    v = bsxfun(@plus,v,curU);
    traceCZ{e} = v;
end
% save figures
preFix = 'eigenVectors_ZSCORE_in_measureSpace';
CON = {'peakPosition','peakValue','Width','asym'};
pairWisePlots(mGn,2,CON,{traceCZ},figure,[],toSave,preFix,toWait);
close all
%% save co-plot
% FIG #4
close all
preFix = 'co_plot_eigenVectors_normalAndZSCORE_in_measureSpace';
CON = {'peakPosition','peakValue','Width','asym'};
pairWisePlots(mGn,2,CON,{traceCZ,traceC},figure,[],toSave,preFix,toWait);
close all
%% 2D plots for eigen curves in u-space
preFix = 'eigenCurves_in_parameterSpace';
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(G,2,LABEL,{ntraceC},figure,[],toSave,preFix,toWait);
close all;
%% 2D plots for eigen curves in u-space - zscore
preFix = 'eigenCurves_in_parameterSpace_zscore';
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(G,2,LABEL,{ztraceC},figure,[],toSave,preFix,toWait);
close all;
%% 2D scatter plots
preFix = 'eigenVectors_in_measureSpace';
CON = {'peakPosition','peakValue','Width','asym'};
pairWisePlots(mGn,2,CON,{ztraceC},figure,[],toSave,preFix,toWait);
%% 3D plots with eigen-curves in u-space
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(G,3,LABEL,{traceC});
CON = {'peakPosition',' peak','Width','asym'};
pairWisePlots(mG,3,CON);
%% trace the unit-vecs in u-space
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(G,2,LABEL,{traceI});
CON = {'peakPosition',' peak','Width','asym'};
pairWisePlots(mG,2,CON);
%% trace both eigen and units
LABEL = {'x0','vf','k','n','len'};
pairWisePlots(G,2,LABEL,{traceI,traceC});
%% scan curves in parameter space to see what they "do" to the regr
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIG #5
close all

preFix = 'parameter_sweeps_along_eigenCurves';
traceToSweep = ntraceC;

preFix = 'parameter_sweeps_along_eigenCurves_Zscore';
traceToSweep = ztraceC;

X = linspace(0,2500,1000);

func = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);
CL = {'r','g','b'};
% conversion factors
conFPS = 30^-1;
conSPH = 60*60;
conPPM = 1463^-1;




% for each curve
for e = 1:min(numel(traceToSweep),3)
    figure
    % for each curve
    for n = 1:2:size(traceToSweep{e},1)
        % eval for velocity
        Y = func(X,traceToSweep{e}(n,:));
        % gradient for regr
        gY = gradient(Y);
        % convert vel
        Y = Y * conFPS * conSPH *conPPM;
        % convert regr
        gY = gY * conFPS * conSPH * 100;
        % co-plot
        yyaxis left
        plot(X,Y,[CL{e} '-']);
        yyaxis right
        plot(X,gY,['k' '-']);
        hold on
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % plot the first
    %%%%%%%%%%%%%%%%%%%%%%%
    n = 1;
    % eval the vel and regr
    Y = func(X,traceToSweep{e}(n,:));
    gY = gradient(Y);
    % convert vel
    Y = Y * conFPS * conSPH *conPPM;
    % convert regr
    gY = gY * conFPS * conSPH * 100;
    % co-plot
    yyaxis right
    plot(X,gY,['m' '-']);
    ylabel('velocity');
    yyaxis left
    plot(X,Y,['m' '-']);
    ylabel('regr');
    xlabel('along axis');
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % plot the last
    %%%%%%%%%%%%%%%%%%%%%%%
    % get the last point on the curve
    n = size(traceToSweep{e},1);
    % eval the vel and regr
    Y = func(X,traceToSweep{e}(n,:));
    gY = gradient(Y);
    % convert vel
    Y = Y * conFPS * conSPH *conPPM;
    % convert regr
    gY = gY * conFPS * conSPH * 100;
     % co-plot
    yyaxis right
    plot(X,gY,['c' '-']);
    yyaxis left
    plot(X,Y,['c' '-']);
    %%%%%%%%%%%%%%%%%%%%%%%
    
   
    fileName = [toSave preFix '_PCAnum_' num2str(e) '.tif'];

    saveas(gca,fileName)
    waitforbuttonpress
    close all
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the axis for the orgSpace
extendPer = 0;
numPoints = 25;
for e = 1:size(G,2)
    minValue = min(G(:,e));
    maxValue = max(G(:,e));
    minValue = minValue*(1-extendPer);
    maxValue = maxValue*(1+extendPer);
    gX{e} = linspace(minValue,maxValue,numPoints);
    dX(e) = (maxValue - minValue)/numPoints;
end
clear gGrid
[gGrid(:,:,:,:,1),gGrid(:,:,:,:,2),gGrid(:,:,:,:,3),gGrid(:,:,:,:,4)] = ndgrid(gX{:});
sz_gGrid = size(gGrid);
gGrid = reshape(gGrid,[prod(sz_gGrid(1:size(G,2))) sz_gGrid(end)]);
zGrid = org_to_z(gGrid);
tic;[mGn] = p2mn(gGrid,3);toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fit a grid to the data - 2D for surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,surC,surU,surE,surL,surERR,surLAM] = PCA_FIT_FULL(mG,2);
[curZ,curU,curSIG] = zscore(mG);
s = diag(surLAM).^.5;
mag = 1.3;
np = 11;
[q1,q2] = ndgrid(linspace(-mag*s(1),mag*s(1),np),linspace(-mag*s(2),mag*s(2),np));
%[q1,q2] = ndgrid(linspace(-mag,mag,np),linspace(-mag,mag,np));
Q = [q1(:) q2(:)];
Q = PCA_BKPROJ(Q,surE,surU);
[responseSurface] = locatePoints(Q,surU,curSIG,zGrid,mGn,z_to_org);
vResponseSurface = responseSurface;
responseSurface = reshape(responseSurface,[np np size(responseSurface,2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fit a grid to the data - 3D for surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toUse = [1 2 4];
[~,surC,surU,surE,surL,surERR,surLAM] = PCA_FIT_FULL(mG,4);
%surU = surU(toUse);
surLAM = diag(surLAM);
surLAM = diag(surLAM(toUse));
surE = surE(:,toUse);
[curZ,curU,curSIG] = zscore(mG);
s = diag(surLAM).^.5;
mag = [20 20 20];
%mag = [1 1 1];
%mag = [2.5 2.5 2.5];
np = [15 15 15];
%np = [3 3 3];
[q1,q2,q3] = ndgrid(linspace(-mag(1)*s(1),mag(1)*s(1),np(1)),linspace(-mag(2)*s(2),mag(2)*s(2),np(2)),linspace(-mag(3)*s(3),mag(3)*s(3),np(3)));
%[q1,q2] = ndgrid(linspace(-mag,mag,np),linspace(-mag,mag,np));
Q = [q1(:) q2(:) q3(:)];
Q = PCA_BKPROJ(Q,surE,surU);
[responseSurface] = locatePoints(Q,surU,curSIG,zGrid,mGn,z_to_org);
vResponseSurface = responseSurface;
responseSurface = reshape(responseSurface,[np size(responseSurface,2)]);
%% plot 3 - zversion for 3D surface
fileName = [toSave 'volumeFly.avi'];
mov1Frames = [toSave 'mov2' filesep];
mkdir(mov1Frames);
mov1 = VideoWriter(fileName);
open(mov1)


dispDims = [1 2 3];
z_G = org_to_z(G);
close all
h = figure;
plot3(z_G(:,dispDims(1)),z_G(:,dispDims(2)),z_G(:,dispDims(3)),'k.');
hold on
for z = 1:size(responseSurface,3)
    for i = 1:size(responseSurface,2)
        tmp = squeeze(responseSurface(:,i,z,:));
        z_tmp = org_to_z(tmp);
        plot3(z_tmp(:,dispDims(1)),z_tmp(:,dispDims(2)),z_tmp(:,dispDims(3)),'r','LineWidth',2);
    end
end

responseSurface = permute(responseSurface,[2 1 3 4]);
for z = 1:size(responseSurface,3)
    for i = 1:size(responseSurface,2)
        tmp = squeeze(responseSurface(:,i,z,:));
        z_tmp = org_to_z(tmp);
        plot3(z_tmp(:,dispDims(1)),z_tmp(:,dispDims(2)),z_tmp(:,dispDims(3)),'g','LineWidth',2);
    end
end
responseSurface = permute(responseSurface,[2 1 3 4]);

responseSurface = permute(responseSurface,[3 2 1 4]);
for z = 1:size(responseSurface,3)
    for i = 1:size(responseSurface,2)
        tmp = squeeze(responseSurface(:,i,z,:));
        z_tmp = org_to_z(tmp);
        plot3(z_tmp(:,dispDims(1)),z_tmp(:,dispDims(2)),z_tmp(:,dispDims(3)),'b','LineWidth',2);
    end
end
responseSurface = permute(responseSurface,[3 2 1 4]);


LABEL = {'x0','vf','k','n','len'};

xlabel(LABEL{dispDims(1)})
ylabel(LABEL{dispDims(2)})
zlabel(LABEL{dispDims(3)})

cameratoolbar('setmode','orbit');
angle = linspace(-90,90,300);
for a = 1:numel(angle)
    view([-angle(a) angle(a)])
    drawnow
    fileName = [mov1Frames num2str(a) '.tif'];
    saveas(gca,fileName);
    I = imread(fileName);

    writeVideo(mov1,I)
end
close(mov1)


%pairWisePlots(zG,3,LABEL,{traceC},h,[1 2 3]);
%% look at the jacobian
zvResponseSurface = org_to_z(vResponseSurface);
alpha = .1;
vecS = [];
for e = 1:size(zvResponseSurface,1)
    % extract the map at each point on the response surface
    [f,jT] = extractMetricsAtP(toMetricSpace_ZT,zvResponseSurface(e,:),.0001*ones(1,size(zvResponseSurface,2)));
    % metric to curvi-linear
    m2c = jT';
    % invert to create the curvilinear to metric
    c2m = inv(m2c);
    
    
    % map eigen vectors in metric space using a different map at each P
    tmp = (c2m*curE(:,1:2));
    
    % get the null space
    tmp = (null(tmp')');
    
    vecS1(e,:) = tmp(1,:);
    vecS2(e,:) = tmp(2,:);
end
vecS1 = zVec_to_org(vecS1);
vecS2 = zVec_to_org(vecS2);
%% find the data on the curviplane
tmpC = PCA_BKPROJ(surC(:,1:2),surE(:,1:2),curU);
[curviData] = locatePoints(tmpC,curU,curSIG,zGrid,mGn,z_to_org);
%% plot 3 - zversion for 2D surface
fileName = [toSave 'surfaceFly.avi'];
mov1Frames = [toSave 'mov1' filesep];
mkdir(mov1Frames);
mov1 = VideoWriter(fileName);
open(mov1)
quivPlot = false;
dispDims = [1 2 3];
qmag = 1;
z_curviData = org_to_z(curviData);
z_G = org_to_z(G);
z_vecS1 = orgVec_to_z(vecS1);
z_vecS2 = orgVec_to_z(vecS2);
close all
h = figure;
plot3(z_G(:,dispDims(1)),z_G(:,dispDims(2)),z_G(:,dispDims(3)),'k.');
hold on
%plot3(z_curviData(:,dispDims(1)),z_curviData(:,dispDims(2)),z_curviData(:,dispDims(3)),'b.');
if quivPlot
    quiver3(zvResponseSurface(:,dispDims(1)),zvResponseSurface(:,dispDims(2)),zvResponseSurface(:,dispDims(3)),z_vecS1(:,dispDims(1)),z_vecS1(:,dispDims(2)),z_vecS1(:,dispDims(3)),qmag,'m')
    %quiver3(zvResponseSurface(:,dispDims(1)),zvResponseSurface(:,dispDims(2)),zvResponseSurface(:,dispDims(3)),z_vecS2(:,dispDims(1)),z_vecS2(:,dispDims(2)),z_vecS2(:,dispDims(3)),2,'c')
end
for i = 1:size(responseSurface,2)
    tmp = squeeze(responseSurface(:,i,:));
    z_tmp = org_to_z(tmp);
    plot3(z_tmp(:,dispDims(1)),z_tmp(:,dispDims(2)),z_tmp(:,dispDims(3)),'r','LineWidth',2);
end
responseSurface = permute(responseSurface,[2 1 3]);
for i = 1:size(responseSurface,2)
    tmp = squeeze(responseSurface(:,i,:));
    z_tmp = org_to_z(tmp);
     plot3(z_tmp(:,dispDims(1)),z_tmp(:,dispDims(2)),z_tmp(:,dispDims(3)),'g','LineWidth',2);
end
responseSurface = permute(responseSurface,[2 1 3]);
LABEL = {'x0','vf','k','n','len'};
xlabel(LABEL{dispDims(1)})
ylabel(LABEL{dispDims(2)})
zlabel(LABEL{dispDims(3)})
%pairWisePlots(zG,3,LABEL,{traceC},h,[1 2 3]);
cameratoolbar('setmode','orbit');
angle = linspace(-90,90,300);
for a = 1:numel(angle)
    view([-angle(a) angle(a)])
    drawnow
    fileName = [mov1Frames num2str(a) '.tif'];
    saveas(gca,fileName);
    I = imread(fileName);

    writeVideo(mov1,I)
end
close(mov1)
%%
angle = linspace(-90,90,500);
for a = 1:numel(angle)
    view([0 angle(a)])
    drawnow
end
angle = linspace(-90,90,500);
for a = 1:numel(angle)
    view([-angle(a) angle(end)])
    drawnow
end
%%

   
splitSpace = [];
for e = 1:size(G,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % map the orginal from G to z
    tmp = org_to_z(G(e,:));
    stmp = org_to_z(curviData(e,:));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % map to measure space
    %tmp = p2mn(tmp);
    %stmp = p2mn(stmp);
    % extract the map at each point on the response surface
    [f_s,jT_s] = extractMetricsAtP(toMetricSpace_ZT,stmp,.0001*ones(1,size(stmp,2)));
    %[f,jT] = extractMetricsAtP(toMetricSpace_ZT,tmp,.0001*ones(1,size(stmp,2)));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parameter to curvi-linear
    p2m = jT_s';
    % invert to create the curvilinear to metric
    m2p = inv(p2m);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % map eigen vectors in metric space using a different map at each P
    eigenVec_inPAtQ = (m2p*curE(:,1:2));
    for v = 1:size(eigenVec_inPAtQ,2)
        eigenVec_inPAtQ(:,v) = eigenVec_inPAtQ(:,v) / norm(eigenVec_inPAtQ(:,v));
    end
    % get the null space
    NULL_inPAtQ = null(eigenVec_inPAtQ');
    
    
    a = eigenVec_inPAtQ'*deltaAtP;
    b = NULL_inPAtQ'*deltaAtP;
    
    
    zztop = zVec_to_org( [a;b]');
    splitSpace(e,:) = zztop;
    mean(splitSpace,1)
    
end
%% look at subtraction
regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;
X = linspace(0,1500,1000);
close all
disp = true;
for e = 1:size(G,1)
    Y = regrFunc(X,G(e,:));
    sY = regrFunc(X,curviData(e,:));
    
    tmp = p2mn(G(e,:));
    stmp = p2mn(curviData(e,:));
    del(e) = norm(curviData(e,:) - G(e,:));
    
    
    
    tmp_u(e) = tmp(2);
    stmp_u(e) = stmp(2);
    
    if disp
        plot(X,Y,'k')
        hold on
        plot(X,sY,'b')
        waitforbuttonpress

        axis([0 1500 0 5*10^-3]);
        drawnow
        pause(.01);
        hold off
    end
    e
end
close all
figure;
plot(del,tmp_u,'.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bubbled to the top from the cell size sand box

found = [];
for tr = 1:100%size(G,1)
    
    
    vel = velFunc(X,G(tr,:));
    regr = regrFunc(X,G(tr,:));
    
    initP = [initPos initLength];
    initP = genCellLength(velFunc,regrFunc,X,G(tr,:),initP,dt,TAU,2.^-cellD);
    cellLength = interp1(initP(:,1),initP(:,2),X,'spline');
    % use the parameter for the velocity but
    tmp_cellLengthPara = G(tr,:);
    % reset the max to be the cell length
    tmp_cellLengthPara(2) = cellLength(end);
    icellGuess = cellLengthFunc(tmp_cellLengthPara);
    plot(X,icellGuess,'g.')
    hold on
    plot(X,cellLength,'r');
    funcy = @(P)norm(cellLengthFunc(P)-cellLength);
    tmp_cellLengthPara = fminsearch(funcy,tmp_cellLengthPara,ops);
    icellGuess = cellLengthFunc(tmp_cellLengthPara);
    plot(X,icellGuess,'k-')
    
    dQ = gradient(vel.*(cellLengthFunc(tmp_cellLengthPara)).^-1,dx);
    figure;
    plot(X*conPPM,dQ*((conPPM)*(conSPH*conFPS)^-1)^-1);
    
    % 
    negCost = 0;
    targetLevel = 10*(conPPM)*(conSPH*conFPS)^-1;
    targetLevel = targetLevel * 100^-1;
    targetLevel = targetLevel;
    
    targetCut = .3*conPPM^-1;
    targetProduction = targetLevel*double(X < targetCut);
    
    targetProduction = imfilter(targetProduction,ones(1,targetKer)/targetKer,'replicate');
    targetProduction = imfilter(targetProduction,fspecial('gaussian',[1,targetKer],targetSig),'replicate');
    ops = optimset('Display','iter','TolX',10^-10,'TolFun',10^-10,'MaxFunEvals',3000,'MaxIter',3000);
    %funcToFix = @(P)norm((cellLengthFunc(P).*gradient(vel.*(cellLengthFunc(P)).^-1,dx) - targetProduction));
    funcToFix = @(P)norm((gradient(vel.*(cellLengthFunc(P)).^-1,dx) - targetProduction));
    
    
    
    found(tr,:) = fminsearch(funcToFix,tmp_cellLengthPara,ops);
    
    found(tr,2) = tmp_cellLengthPara(2);
    
    cellLengthFunc = @(P)velFunc(X,P);
    
    Q = cellLengthFunc(found(tr,:));
    close all
    plot(X*conPPM,Q*conPPM)
    
    dQ = gradient(vel.*(cellLengthFunc(found(tr,:))).^-1,dx);
    figure;
    plot(X*conPPM,dQ*((conPPM)*(conSPH*conFPS)^-1)^-1);
    hold on
    plot(X*conPPM,targetProduction*((conPPM)*(conSPH*conFPS)^-1)^-1,'r')
    waitforbuttonpress
      
    %initP = [initPos initLength];
    %what = genCellLength(velFunc,regrFunc,X,G(tr,:),initP,dt,TAU,2.^-abs(dQ));
end
%%

kp = find(~isoutlier(found(:,end)));
for e = 1:size(G,2)
    figure
    plot(G(kp,e),found(kp,2)/1463,'k.')
    waitforbuttonpress
end
%%
for e = 1:numel(fidx)
    cL = cellLengthFunc(found(kp(e),:));
    plot(cL);
    waitforbuttonpress
end
%%
close all
tr = 1;
    
X = linspace(0,2500,10000);
dx = mean(diff(X,1,2));
initPos = 10*1463^-1;
initPos = 10;
initLength = .20;
initLength = .020*1463;

mag = 1;
initP = [initPos initLength];
dt = 1;
TAU = 24*5;


conFPS = 30^-1;
conSPH = 60*60;

dt = dt*conSPH*conFPS;
TAU = TAU*conSPH*conFPS;



velFunc = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);
regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;



vel = velFunc(X,G(tr,:));
regr = regrFunc(X,G(tr,:));

[cellD] = genCellD(X,[500 200 .01]);



%cellD = zeros(size(X));
figure;
plot(X,cellD);
waitforbuttonpress
close all

%{
Y = double(X < para(1));
%Y = imfilter(Y,ones(1,para(2))/para(2),'replicate');
%Y = imfilter(Y,ones(1,para(2))/para(2));
Y = imfilter(Y,fspecial('gaussian',[1 2*para(2)],para(2)),'replicate');
Y = para(3) * Y * max(Y(:))^-1;
%}


cellD = zeros(size(cellD));

%%
%%%%%%%%%%%%%%%%%%%%
close all
cellLengthFunc = @(P)velFunc(X,P);
delta = @(P)norm(cellLengthFunc(P) - cellLength);
initP = G(tr,:);
initP(2) = cellLength(end);
tmp_cellLengthPara = fminsearch(delta,initP);
Y = cellLengthFunc(tmp_cellLengthPara);
plot(X,Y,'g.');
hold all
plot(X,cellLength);

%%
%f = fit(X',cellDensity','exp1');
%a = f.a;
%b = f.b;
c = 1;
cellDensity = a*exp(X*b*c);
figure;
plot(X,cellDensity.^-1);
%%
close all
vel = velFunc(X,G(tr,:));
regr = regrFunc(X,G(tr,:));
per = linspace(1,.75,10);


for e = 1:20
    
    
    
    
    
    
    
    initP = [initPos initLength];
    initP = genCellLength(velFunc,regrFunc,X,G(tr,:),initP,dt,TAU,2.^-cellD);
    cellLength = interp1(initP(:,1),initP(:,2),X,'spline');
    
    tmp_cellLengthPara = G(tr,:);
    tmp_cellLengthPara(2) = cellLength(end);
    
   
    
    negCost = 0;
    targetLevel = 10^-5;
    targetCut = 300;
    targetKer = 901;
    targetSig = 300;
    targetLevel = 10*(conPPM)*(conSPH*conFPS)^-1;
    targetCut = .3*conPPM^-1;
    targetProduction = targetLevel*double(X < targetCut);
    targetProduction = imfilter(targetProduction,ones(1,targetKer)/targetKer,'replicate');
    targetProduction = imfilter(targetProduction,fspecial('gaussian',[1,targetKer],targetSig),'replicate');
    %targetProduction = imfilter(targetProduction,fspecial('gaussian',[1,targetKer],targetSig));
    ops = optimset('Display','iter','TolX',10^-10,'TolFun',10^-10,'MaxFunEvals',3000,'MaxIter',3000);
    funcToFix = @(P)norm(sum(negCost*((gradient(vel.*(cellLengthFunc(P)).^-1,dx))<0))...
        +(gradient(vel.*(cellLengthFunc(P)).^-1,dx) - targetProduction));
    found = fminsearch(funcToFix,tmp_cellLengthPara,ops);
    
    
   
    Q = cellLengthFunc(found);
    plot(X,Q)
    dQ = gradient(vel.*(cellLengthFunc(found)).^-1,dx);
    figure;
    plot(X,dQ);
    hold on
    plot(X,targetProduction,'r')
    waitforbuttonpress
    
    
    convect = vel.*gradient(cellDensity,dx);
    stretch = regr.*cellDensity;
    cellProduction = gradient(vel.*cellDensity,dx);
    plot(X,convect,'r');hold on;
    plot(X,stretch,'g');
    plot(X,cellProduction,'b');
    waitforbuttonpress
    hold off
    pause(.1)
    
    
    
    %{
    waitforbuttonpress

    figure;
    yyaxis left
   
    
    plot(convect,'r');
    hold on
    plot(stretch,'k')
    %legend({'Convective','Stretch'});
    %}
    %{
    ylabel('Convective');
    yyaxis right
    plot(convect + stretch,'k');
    hold on
    plot(cellProduction,'g--');
    ylabel('Stretch');
    hold off
    drawnow
    pause(.1)
    %}
    %waitforbuttonpress
end
    

%%
for e = 1:50
    
    initP = [initPos initLength];
    
    [initP] = genCellLength(velFunc,regrFunc,X,G(tr,:),initP,dt,TAU,2.^-cellD);
    cellLength = interp1(initP(:,1),initP(:,2),X,'spline');
    rm = any(isnan(initP),2);
    initP(rm,:) = [];
    %xlabel('days')
    %ylabel('cell size')
    
    
    
    
    figure;
    yyaxis left
    plot(X,cellLength);
    ylabel('Cell Length');
    yyaxis right
    plot(X,regr);
    ylabel('REGR');
    cellDensity = cellLength.^-1;
    cellProduction = gradient(vel.*cellDensity,dx);
    %waitforbuttonpress
    
    
    cellD = cellProduction;
    
    figure;
    yyaxis left
    plot(X,2.^-cellD);
    ylabel('Cell Production');
    yyaxis right
    plot(X,regr);
    ylabel('REGR');
    %waitforbuttonpress
    
    
    figure;
    yyaxis left
    convect = vel.*gradient(cellDensity,dx);
    stretch = regr.*cellDensity;
    plot(convect,'r');
    hold on
    plot(stretch,'k')
    %legend({'Convective','Stretch'});
   
    ylabel('Convective');
    yyaxis right
    plot(convect + stretch,'k');
    hold on
    plot(cellProduction,'g--');
    ylabel('Stretch');
    waitforbuttonpress
    
    drawnow
    close all
end



%%
figure;
plot(2.^-cellProduction);
%% JUNK YARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;
tr = 6;


paraP = uean_AtPara;
pw = repmat(paraP,[N+1 1])';

[X,W] = ndgrid(linspace(0,2500,2500),linspace(-N/2,N/2,N+1));
for w = 1:size(pw,2)
    VEL_X(:,w) = velFunc(X(:,1),pw(:,w));
    VEL_W(:,w) = zeros(size(VEL_X(:,w)));
end
mesh(VEL_X)

%%

traceVecSet = zcurE(:,1);
nTrace = (size(W,2)-1)/2;
maxV = .5;
alpha = maxV/nTrace;
for e = 1%1:size(traceVecSet,2)
    dv = traceVecSet(:,e);
    traceInLat1_pos = org_to_z(uean_AtPara);
    [traceInLat1_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_pos,dv,nTrace,alpha);
    traceInLat1_neg = org_to_z(uean_AtPara);
    [traceInLat1_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_neg,-dv,nTrace,alpha);
    traceInLat1_neg(1,:) = [];
    delta = z_to_org([flip(traceInLat1_pos,1);traceInLat1_neg]);
end
delta = delta';
%%
close all
[~,angV] = curl(VEL_W,VEL_X);
mag1 = .5;
%lam = linspace(-mag1,mag1,size(W,2));
%delta = curE(:,1)*lam;
pw_p = pw + delta;
for w = 1:size(pw,2)
    VEL_X_p(:,w) = velFunc(X(:,1),pw_p(:,w));
    VEL_W_p(:,w) = zeros(size(VEL_X(:,w)));
end
mesh(VEL_X_p)
[~,angV] = curl(VEL_W,VEL_X);
%% JUNK YARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% test bed for zoom
W.numP = [10 10 10 10];
W.width = .1*std(G,1,1);
zoomP = mean(G,1);
zoomOnPoint(@(p)p2mn(p),zoomP,W);
%% test bed for point finding
clear p
targetQ = curU;
disToP = @(p)sum(bsxfun(@minus,mGn,p).^2,2);
isoNewP = @(d)find(d(:) == min(d(:)));
disToTarget = @(X)sum(bsxfun(@minus,p2mn(X,11),targetQ).^2,2);
d = disToP(curU);
didx = isoNewP(d);
initP = gGrid(didx,:);
globalNumP = 5*ones(1,size(gGrid,2));
floatLevel = 4;
loopMax = 20;

p = initP;
F = disToTarget;
W.numP = globalNumP;
W.width = dX;
L.totalLoop = 0;
L.resolutionValue = 0;
L.loopCompression = .1;
L.updateCurrentResolution = @(p,history)-log10(max(abs(history - p)));
L.history = [];
L.halt = @(loop,res)res >= floatLevel | loop > loopMax;
isoNewP = isoNewP;
[p,retL] = ifunc(p,F,W,isoNewP,L);
q = p2mn(q);
targetQ
q
%}




%{
%% this might not be needed - the idea was to trace out the grid
% instead i am trying the below method
gridTrace(uean_AtPara,toMetricSpace_ZT,org_to_z,z_to_org,curE(:,1),curE(:,2));
%}