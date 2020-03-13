%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% locate on network drives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% locate all image files
[FileList] = locateImages();
lFileList = lower(FileList);
%% locate only files that are named with numbers
[nFileList] = locateNumberedImageSets();
%% locate by person
personList = lower({'Ashley'});
personList = lower({'nate'});
pIDX = contains(lFileList,personList);
pFileList = FileList(pIDX);
numel(pFileList)
%% scan for locations
nFileList = fastOrderintoSets(nFileList);
%% report on nFileList
nSZ = [];
fIndex = {};
for e = 1:numel(nFileList)
    nSZ(e) = numel(nFileList{e});
    fIndex{e} = nFileList{e}{1};
end
nSZ
numel(nFileList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seed and grain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE: SAG-Seed and Grain Algorithm Layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make init vars layer
initContext_SAG = initLayer('SAG: Initialization Layer');
% make the seed image read layer
seedReadLayer = readLayer('generalSeedReaderLayer');
seedReadLayer.isConfigured = true;
seedReadLayer.setType('bvRGBsimpleImageLoader');
rgb2hsvLayer = colorConvertLayer('RGB->HSV','rgb>hsv');
cc = functionHandleLayer('Custom Seed and Grain Layer',...
    @(X,Y,Z,A)seedAndGrain_customCode_block1(X,Y,Z,A));
seedLayerStack = [initContext_SAG;seedReadLayer;rgb2hsvLayer;cc];
SAG_algorithm = algorithm('SAG-Seed and Grain',seedLayerStack);
SAG_algorithm.setAuthor('Nathan Miller');
SAG_algorithm.setVersion('1');
SAG_algorithm.setDescription(['Method measures seed and grain traits. From images taken with camera or ' ...
    'flatbed scanner.  Method can be used wth light seeds on dark background or dark seeds '...
    'on light background.']);
%% search for terms in the file name for seed images
containKeyWords = {'seed'};
notContainKeyWords = {'logo','results','media','figure','root','seedsource','seedling','swell','swell','kernel','corn','maize','hypo','return'};
seedIDX = find(contains(lFileList,containKeyWords) & ~contains(lFileList,notContainKeyWords));
seedFileList = FileList(seedIDX);
%% create solid test set
containKeyWords = {'seed'};
notContainKeyWords = {'phytomorph','logo','results','media','figure','root','seedsource','seedling','swell','swell','kernel','corn','maize','hypo','return'};
seedIDX = find(contains(lFileList,containKeyWords) & ~contains(lFileList,notContainKeyWords));
seedFileList_smallTest = FileList(seedIDX);
%% create random test set
N = 100;
rFileList = seedFileList(randperm(numel(seedFileList)));
rFileList = rFileList(1:N);
%% add images to random test procedure
[targetList] = networkTester.addImageToTestProcedure(rFileList,SAG_algorithm,'randomSet','link');
%% create standard test set
N = 10;
sdFileList = seedFileList(randperm(numel(seedFileList)));
sdFileList = sdFileList(1:N);
%% add images to random test procedure
[targetList] = networkTester.addImageToTestProcedure(sdFileList,SAG_algorithm,'standardSet','link');
%% list images in random test procedure
[tFileList] = networkTester.listImagesForTestProcedure(SAG_algorithm,'standardSet');
%% set standards for standard test procedure
networkTester.setStandard(SAG_algorithm,'standardSet');
%% run test
networkTester.test(SAG_algorithm,'standardSet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seed and grain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MECKA - cob data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make Maize Cob Layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image resolution
scanResolution = '1200';
% set to default value of 1200
defaultResolution = 1200;
scanResolution = StoN(scanResolution);
% compute proportion of resolution over default
fracDpi = scanResolution/defaultResolution;
% set to default value of .25
checkBlue_scaleFactor = .25; 
% set to default value of 100
addcut = 150;
addcut = round(addcut*fracDpi);
% set to default value of 600
baselineBlue = 600;
baselineBlue = round(baselineBlue*fracDpi);
% set to default value of 10^6
defaultAreaPix = 10^6;
defaultAreaPix = round(defaultAreaPix*fracDpi);
% set to default value of 300
rho = 300;
rho = round(rho*fracDpi);
% set to default value of 70
colRange1 = 70;
% set to default value of 166
colRange2 = 166;
% set to default value of 50
fill = 50;
% set number of objects
noe = 3;
rPath = [];
rawImage_scaleFactor = 1;
toSave = false;
toDisplay = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE: MECKA-Cob layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make init vars layer
initContext_cMECKA = initLayer('MECKA:Cob:Initialization Layer');
MECKA_cobLayer = functionHandleLayer('Custom MECKA cob layer.',...
    @(X,Y)singleCobImage_octerine...
    (X,noe,Y,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,...
    defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay));


MECKAcob_LayerStack = [MECKA_cobLayer];
% make algorithm from layers
MECKAc_algorithm = algorithm('MECKA-Maize Ear,Cob and Kernel-Cob Version',MECKAcob_LayerStack);
MECKAc_algorithm.setAuthor('Nathan Miller');
MECKAc_algorithm.setVersion('1');
MECKAc_algorithm.setDescription(['Method measures maize cobs including width profile, length, and color.']);

%% search for terms in the file name for seed images
containKeyWords = {'cobdata'};
notContainKeyWords = {'what'};
cobIDX = find(contains(lFileList,containKeyWords) & ~contains(lFileList,notContainKeyWords));
cobFileList = FileList(cobIDX);
%% create standard test set
N = 10;
sdFileList = cobFileList(randperm(numel(cobFileList)));
sdFileList = sdFileList(1:N);
%% add images to random test procedure
[targetList] = networkTester.addImageToTestProcedure(sdFileList,MECKAc_algorithm,'standardSet','link');
%% list images in random test procedure
[tFileList] = networkTester.listImagesForTestProcedure(MECKAc_algorithm,'standardSet');
%% set standards for standard test procedure
networkTester.setStandard(MECKAc_algorithm,'standardSet');
%% run test
networkTester.test(SAG_algorithm,'standardSet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MECKA - cob data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 

oPath = '/home/nate/10_26_trials/';
fileName = '/home/nate/10_26_trials/GS_test_COL_4800dpi_rgb028.tif';
fileName = '/iplant/home/turnersarahd/Clement_drought_scans/Normal/035255243056';
fileName = '/home/nate/Blue_wheat.jpg';


FilePath = '/home/nate/10_26_trials/Calibration Images/';
FileList = {};
FileExt = {'tiff','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
   

%% look for hypocotyal images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hypoIDX = contains(fIndex,{'hypo'},'IgnoreCase',true) & ~contains(fIndex,{'julian','matlab'},'IgnoreCase',true);
%hypoIDX = contains(fIndex,{'dark'},'IgnoreCase',true) & ~contains(fIndex,{'julian','matlab'},'IgnoreCase',true);
%hypoIDX = contains(fIndex,{'wt'},'IgnoreCase',true) & ~contains(fIndex,{'kinematics','root','julian','matlab'},'IgnoreCase',true);
hFileList = fIndex(hypoIDX);
hFileList
%%
for e = 1:numel(hFileList)
    I = imread(hFileList{e});
    imshow(I,[]);
    drawnow
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mass download for seed size tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target = '/mnt/spaldingdata/nate/massMaize/';
source = '/iplant/home/nhaase/%';
imassDownload(source,target,'tif',12)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mass download for seed size tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target = '/mnt/spaldingdata/nate/Craines/';
source = '/iplant/home/ssslab/ScanImages/%';
imassDownload(source,target,'tif',12)
%% REVIST THE CREATING OF THE TEST PROCEDURE IMAEGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:numel(seedIDX)
    seedFileList{e} = networkTester.addImageToTestProcedure(FileList{seedIDX(e)},'SAG','randomSet','link');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select kernel list(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
containKeyWords = {'kernel'};
notContainKeyWords = {'logo','output','snowbird','results','swell','junk','figure','outport','return'};
kernelIDX = find(contains(lFileList,containKeyWords) & ~contains(lFileList,notContainKeyWords));
kernelFileList = FileList(kernelIDX);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% view the seedFileList
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
for e = 1:numel(kernelFileList)
    I = imread(kernelFileList{e});
    imshow(I,[]);
    drawnow
end
%%
for e = 1:numel(FileList)
    lower(FileList)
end
w = find(contains(FileList,'seed'))

%%
testFunction = 'measureArabidopsisSeed_BK';
commenter = autoCommenter();
[funcHandle] = commenter.moniterVars(testFunction);
global totalMem;
global totalCount;
global lineN;

lineN = [];
totalMem = [];

oPath = '/home/nate/output/';

localM = [];
totalCount = [];
maxLine = [];
topN = 4;

for e = 1:numel(seedList)
    
    funcHandle(seedList{e},oPath);
    
    dM = diff(totalMem);
    [mm,midx] = sort(dM,'descend');
    %[mm,midx] = max(dM);
    
    
    maxLine(e,:) = lineN(midx(1:topN));
    
    
    featureM = [max(totalMem) mean(totalMem) std(totalMem)];
    
    
    
    
    localM = [localM;featureM];
    e
    
    
    totalMem = [];
    totalCount = [];
    lineN = [];
    
    
    plot(localM);
    drawnow
    maxLine
end




%%
global totalMem 
MX = max(log10(totalMem));
GB = 9;
MB = 6;
KB = 3;
B = 1;
close all
LEG = {'Memory Use'};
plot(log10(totalMem))
hold all
if 2*MX > KB;plot(KB*ones(size(totalMem)));LEG{end+1} = 'KB';end
if 2*MX > MB;plot(MB*ones(size(totalMem)));LEG{end+1} = 'MB';end
if 2*MX > GB;plot(GB*ones(size(totalMem)));LEG{end+1} = 'GB';end
%legend(LEG)
axis([0 numel(totalMem) 0 10]);

for e = 2:5
    plot((GB+log10(e))*ones(size(totalMem)));
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make Maize Cob Layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image resolution
scanResolution = '1200';
% set to default value of 1200
defaultResolution = 1200;
scanResolution = StoN(scanResolution);
% compute proportion of resolution over default
fracDpi = scanResolution/defaultResolution;
% set to default value of .25
checkBlue_scaleFactor = .25; 
% set to default value of 100
addcut = 150;
addcut = round(addcut*fracDpi);
% set to default value of 600
baselineBlue = 600;
baselineBlue = round(baselineBlue*fracDpi);
% set to default value of 10^6
defaultAreaPix = 10^6;
defaultAreaPix = round(defaultAreaPix*fracDpi);
% set to default value of 300
rho = 300;
rho = round(rho*fracDpi);
% set to default value of 70
colRange1 = 70;
% set to default value of 166
colRange2 = 166;
% set to default value of 50
fill = 50;

% make init vars layer
initContext_MECKAcob = initLayer('MECKA:Cob Initialization Layer');
%%
rawImage_scaleFactor = 1;
scanResolution = 1200;
toDisplay = 1;
toSave = 1;
remotePath = [];
oPath = '/home/nate/Downloads/';
numberOfObjects = 3;

fileName = '/iplant/home/kmichel/maizeData/cobData/17-06-02/WISN16_110_imageData_cob.tif';
fracDpi = 1;
% set to default value of 10^6
defaultAreaPix = 10^6;
defaultAreaPix = round(defaultAreaPix*fracDpi);
% set to default value of 300
rho = 300;
rho = round(rho*fracDpi);
% set to default value of 70
colRange1 = 70;
% set to default value of 166
colRange2 = 166;
% set to default value of 50
fill = 50;

 % set to default value of 1200
    defaultResolution = 1200;
    scanResolution = StoN(scanResolution);
    % compute proportion of resolution over default
    fracDpi = scanResolution/defaultResolution;
    % set to default value of .25
    checkBlue_scaleFactor = .25; 
    % set to default value of 100
    addcut = 150;
    addcut = round(addcut*fracDpi);
    % set to default value of 600
    baselineBlue = 600;
    baselineBlue = round(baselineBlue*fracDpi);

    
    
    
singleCobImage(fileName,numberOfObjects,oPath,remotePath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay);
           

%%
oldPath = path;
%%

%%
testFunction = 'singleCobImage_octerine';
commenter = autoCommenter();
[newFunctionName,newFunctionFile,funcHandle] = commenter.moniterVars(testFunction);
%%
clear global
funcHandle(fileName,numberOfObjects,oPath,remotePath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay);
%%
global totalMem 
MX = max(log10(totalMem));
GB = 9;
MB = 6;
KB = 3;
B = 1;
close all
LEG = {'Memory Use'};
plot(log10(totalMem))
hold all
if 2*MX > KB;plot(KB*ones(size(totalMem)));LEG{end+1} = 'KB';end
if 2*MX > MB;plot(MB*ones(size(totalMem)));LEG{end+1} = 'MB';end
if 2*MX > GB;plot(GB*ones(size(totalMem)));LEG{end+1} = 'GB';end
legend(LEG)
axis([0 numel(totalMem) 0 10]);

%%
convertFuntionForm('singleCobImage_octerine_TEST')
%%
%%
dataPoolPort.makeDataPool('images/thumbNails');
dataPoolPort.removeDataPool('images/thumbNails/');
%%





%%
%nTester = networkTester();
%nTester.test(SAG_algorithm,'randomSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% view the seedFileList
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
for e = 1:numel(seedFileList)
    I = imread(seedFileList{e});
    imshow(I,[]);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test for the SAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% "manual" test
timingBlock('clear')
FilePath = '/home/nate/10_26_trials/Calibration Images/';
FileList = {};
FileExt = {'tiff','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
fluid = digitalFluid('randTestFluid',{FileList(1),'./placeHere/'});
tic;fluid = fluid.flow(SAG_algorithm);etim = toc;
%% place network in testing device
nTester = networkTester();
nTester.test(SAG_algorithm,'randomSet');