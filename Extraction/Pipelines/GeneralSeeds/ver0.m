FilePath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/ColorAnalysis/';
FileSet = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FileSet = sfdig(FilePath,FileSet,FileExt,1);
%% view the data
close all
skip = [100 100];
X = [];
Y = {};
H = [];
M = [];
cnt = 1;
hbins = linspace(0,1,255);
h1 = figure;
h2 = figure;
disp = false;
darkT = 20;
% for each set
for s = 1:numel(FileSet)
    FileList = FileSet{s};
    
    tX = [];tY = {};tH = [];
    % jfor each image in the set
    parfor e = 1:numel(FileList)
        
        [fld] = extractFieldFromPath(FileList{e},1);
        
        info = imfinfo(FileList{e});
        ROWS = [1 skip(1) info.Height];
        COLS = [1 skip(2) info.Width];
        I = double(imread(FileList{e},'PixelRegion',{ROWS,COLS}))/255;
        
        
        % find the black connected to the edge - the black background
        %%%%%%%%%%%%%%%%%%%%
        LAB = rgb2lab(I);
        dark = LAB(:,:,1) < darkT;
        darkEdge = (imclearborder(dark) == 0) & (dark == 1);
        colorArea = ~darkEdge;
        cidx = find(colorArea);
        
        % build up the histogram for the image only at the
        % color-background-area
        %%%%%%%%%%%%%%%%%%%%
        tmpH = [];
        for k = 1:3
            tmpI = I(:,:,k);
            tmpH(:,k) = hist(tmpI(cidx),hbins);
            tmpH(:,k) = tmpH(:,k) / sum(tmpH(:,k));
        end
        
        % display code
        %%%%%%%%%%%%%%%%%%%%
        if disp
            % show the image
            figure(h1);
            imshow(I,[]);
            title(fld);

            % plot the histogram
            figure(h2);
            plot(tmpH);
            drawnow
        end
        
        % store X=imagedata,Y=imagetype,H=histogram,M=mask
        %%%%%%%%%%%%%%%%%%%%
        tX(:,:,:,e) = I;
        tY{e} = fld;
        tH(:,:,e) = tmpH;
        tM(:,:,e) = colorArea;
    end
    
    s
    
    X = cat(4,X,tX);
    Y = cat(2,Y,tY);
    H = cat(3,H,tH);
    M = cat(3,M,tM);
end
%% look at the backgrounds and histograms for each type
UQ = unique(Y);
close all
CL = {'r','g','b'};
for u = 1:numel(UQ)
    fidx = strcmp(Y,UQ{u});
    tmpH = mean(H(:,:,fidx),3);
    for k = 1:3
        plot(hbins,tmpH(:,k),CL{k}); 
        hold on
    end
    title(UQ{u});
    hold off
    axis([0 1 0 .5])
    waitforbuttonpress
end
%% make average image for each type
UQ = unique(Y);
close all
toMap = [];
CL = {'r','g','b'};
for u = 1:numel(UQ)
    % find the trials of type-u
    fidx = strcmp(Y,UQ{u});
    % average by type
    uI = mean(X(:,:,:,fidx),4);
    % store average
    uIstore(:,:,:,u) = uI;
    % make mask
    samM = all(M(:,:,fidx),3);
    % find mask == 1
    midx = samM(:) == 1;
    % sample color space at mask
    samV = [];
    for k = 1:3
        tI = uI(:,:,k);
        samV = [samV tI(midx)];
    end
    [U,E,L] = PCA_FIT_FULLws(samV,3);
    mapLevel_in{u} = @(x)PCA_REPROJ(x,E,U);
    mapLevel_out{u} = @(x)PCA_BKPROJ(x,E,U);
    toMap(:,:,u) = [E];
    
    probMap{u} = @(x)mvnpdf(x,U,cov(samV));
    maxValue(u) = probMap{u}(U);
    %toMap(:,:,u) = [[U' bsxfun(@plus,E,U')'];[0 0 0 1]];
    %toMap(:,:,u) = [[zeros(size(U')) bsxfun(@plus,E,U')'];[1 1 1 1]];
    % show
    out = flattenMaskOverlay(uI,samM > .8);
    imshow(out,[]);
    title(UQ{u});
    pause(.5);
    %waitforbuttonpress
end
%% map between color spaces
n = nchoosek(1:size(toMap,3),2);
n = [n;flip(n,2)];
map = {};
for e = 1:size(n,1)
    colorMap = toMap(:,:,n(e,1)) \ toMap(:,:,n(e,2));
    testV = abs(toMap(:,:,n(e,1))*colorMap - toMap(:,:,n(e,2))) < .001;
    all(testV(:))
    
    map{n(e,1),n(e,2)} = colorMap;
    
    colorMapper{n(e,1),n(e,2)} = @(x)mapLevel_out{n(e,2)}((colorMap*(mapLevel_in{n(e,1)}(x))')');
    
    
end
for e = 1:4
    map{e,e} = eye(4);
    colorMapper{e,e} = @(x)x;
end
%% try matching a full image using PCA mapper
close all
sourceSet = 1;
targetSet = 3;
n = 5;
darkT = 20;
skip = [10 10];
info = imfinfo(FileList{e});
ROWS = [1 skip(1) info.Height];
COLS = [1 skip(2) info.Width];
I = double(imread(FileSet{sourceSet}{n},'PixelRegion',{ROWS,COLS}))/255;
[colorArea,cidx] = findColorBackgroundArea(I,darkT);
toUseFunc = colorMapper{sourceSet,targetSet};
[aI] = rasterColorImage(I,toUseFunc,cidx);
imshow([I,aI],[]);
%% try to map data to arabidopsis seed method - testing here
close all

out = {};
colorPatch = {};
count = [];
oPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/ColorAnalysis/return/';
for sourceSet = 1:4
    for targetSet = 1:4

        s_fld = extractFieldFromPath(FileSet{sourceSet}{1},1);
        t_fld = extractFieldFromPath(FileSet{targetSet}{1},1);
        backgroundType = {};
        toProcess = false;
        
        parfor n = 1:20
            try
                n
              
                fileName = FileSet{sourceSet}{n};
                [~,nm] = fileparts(fileName);


                toPath = [oPath '{source_' s_fld '}' '{target_' t_fld '}' filesep];
                mmkdir(toPath);


                darkT = 20;
                skip = [5 5];
                info = imfinfo(fileName);
                ROWS = [1 skip(1) info.Height];
                COLS = [1 skip(2) info.Width];

  

                I = double(imread(fileName,'PixelRegion',{ROWS,COLS}))/255;
                %{
                if toProcess

                    networkSkip = [100 100];
                    ROWS = [1 networkSkip(1) info.Height];
                    COLS = [1 networkSkip(2) info.Width];
                    nI = double(imread(FileSet{sourceSet}{n},'PixelRegion',{ROWS,COLS}))/255;
                    backgroundType{n} = net.classify(nI);



                    [colorArea,cidx] = findColorBackgroundArea(I,darkT);
                    toUseFunc = colorMapper{sourceSet,targetSet};
                    [aI] = rasterColorImage(I,toUseFunc,cidx);
                    [pI] = rasterColorImage(aI,probMap{targetSet},cidx,maxValue(targetSet));


                    logSig = -log(pI+eps*10^-100);

                    pIb = bindVec(logSig);
                    pIi = cat(3,pIb,pIb,pIb);

                    [out{sourceSet,targetSet,n},colorPatch{sourceSet,targetSet,n},count(sourceSet,targetSet,n)] = measureGeneralSeeds(I,pIb);

                    oFile = [toPath nm '{predict_' char(backgroundType{n}) '}.jpg'];
                    imwrite(out{sourceSet,targetSet,n},oFile);
                end
                %}
                
                oFile = [toPath nm '{type_original}.jpg'];
                imwrite(I,oFile);
                
            catch ME
                ME
            end
            
        end
        
        
        %{
        ROWS = [1 skip(1) info.Height];
        COLS = [1 skip(2) info.Width];
            
            
        for n = 1:20
            I = double(imread(FileSet{sourceSet}{n},'PixelRegion',{ROWS,COLS}))/255;
            cA = imcrop(I,[1 1 600 600]);
            %cB = imcrop(aI,[1 1 600 600]);
            %cC = imcrop(pIi,[1 1 600 600]);
            cD = imcrop(out{sourceSet,targetSet,n},[1 1 600 600]);
            close all
            imshow([cA,cD],[]);
            title([s_fld '-->' t_fld '-->' char(backgroundType{n})]);
            pause(.5);
        end
        %}
        
        
    end
end
%% try matching a full image using built-in histogram matching
sourceSet = 3;
targetSet = 1;
n = 5;
skip = [10 10];
info = imfinfo(FileList{e});
ROWS = [1 skip(1) info.Height];
COLS = [1 skip(2) info.Width];
I = double(imread(FileSet{sourceSet}{n},'PixelRegion',{ROWS,COLS}))/255;
aI = imhistmatch(I,uIstore(:,:,:,targetSet),255);
imshow([aI I],[]);
%% set up the data for training network on backgroun type - categorical
szX = size(X);
UQ = unique(Y);
c = cvpartition(numel(Y),'HoldOut',.2);
xTrain = X(:,:,:,c.training);
yTrain = categorical(Y(c.training));
xTest = X(:,:,:,c.test);
yTest = categorical(Y(c.test));
%% train the network
layers = [imageInputLayer(szX(1:3),'Normalization','none');
          convolution2dLayer(11,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(numel(UQ));
          softmaxLayer();
          classificationLayer()];
options = trainingOptions('sgdm','MaxEpochs',1000,'InitialLearnRate',0.001,...
        'ExecutionEnvironment','parallel','Plots','training-progress',...
        'ValidationData',{xTest,yTest'});
    
net = trainNetwork(xTrain,yTrain',layers,options);
%% look at registration
% results - does not appear that the images are registered
close all
skip = [10 10];
info = imfinfo(FileSet{1}{1});
ROWS = [1 skip(1) info.Height];
COLS = [1 skip(2) info.Width];
N = 10;
I = double(imread(FileSet{1}{N},'PixelRegion',{ROWS,COLS}))/255;
J = double(imread(FileSet{2}{N},'PixelRegion',{ROWS,COLS}))/255;
K = double(imread(FileSet{3}{N},'PixelRegion',{ROWS,COLS}))/255;
I = rgb2gray(I);
J = rgb2gray(J);
K = rgb2gray(K);
imshow(cat(3,I,J,K),[]);
%%
      