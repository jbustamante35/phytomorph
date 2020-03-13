%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% opti functions look like they are the meta fusion layers
%%%
% This file is for the sequence of functions needed to extract, train and
% pubish a pipeline which detects "dots" in an image
%% load the main program for viewing and edits
websave(['./sorghumNNapp.mat'],'https://de.cyverse.org/dl/d/A2E6EFA2-37C3-4EC6-9335-302C29F9B2CC/metaLayers_sorghum.mat');
obj = load(['./sorghumNNapp.mat']);
%%
inlineFeatureExtractionProgram = ['./generalizeFeatureExtractor.mat'];
% save generalized loader from irods
websave(inlineFeatureExtractionProgram,'https://de.cyverse.org/dl/d/5C41DB00-9619-4FB0-A18B-66757FC1503A/generalizeFeatureExtractor.mat');
fe = load(inlineFeatureExtractionProgram);
%% gather all maizeData
% this will only comb through the directory and find the maize
dataPath = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        maizeFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%% gather all sorghumData
dataPath = ['/iplant/home/leakey_cyverse/sorghumData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        sorFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%%
[readerHistogram_sor] = gatherParametersForReader(sorFileList(1:100),'histogramNormalize',false);
[readerHistogram_maize] = gatherParametersForReader(maizeFileList(1:100),'histogramNormalize',false);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate operational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I = imread(sorFileList{1});
I = imread(maizeFileList{1});
%I = imread(setFileList{1});
imageSZ = size(I);
borderSize = [40 40];
[pointSet,pointSetSize] = genIXgrid(imageSZ,[1 1],borderSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample domain-1 disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opDomain = {};
domainSize = {};
R = [0 40];
N = [41 250];
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
X = n1.*cos(n2+pi);
Y = n1.*sin(n2+pi);
opDomain{1} = [X(:) Y(:)];
domainSize{1} = size(X);
testT{1} = opDomain{1};



freezeTensor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package input(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trimValues = R;
pointSetData = pointSet;
domainData{1} = opDomain;
domainData{2} = domainSize;
domainData{3} = trimValues;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make input list(s)
% make single input
inputPoint = {sorFileList{1},pointSet(1,:),opDomain,domainSize};
for e = 1:10
    inputPoints{e} = {sorFileList{1},pointSet(e,:),opDomain,domainSize};
end
for e = 1:size(pointSet,1)
    TinputPoints{e} = {sorFileList{1},pointSet(e,:),opDomain,domainSize};
end

inputPoints = inputPoints';
TinputPoints = TinputPoints';

%% make network layers
% make the read layer
sorReadLayer = readLayer('sorghumStomataReader');
sorReadLayer.setType('histogram');
% make the sample layer
sampleLayer = samplerLayer('SampleStomata');
% make and configure the fft layer
FFTLayer = fftLayer('sorghumStomataFFT');
FFTLayer.setN(35);
% make and configure the compression layer
FFTCompressionLayer = compressionLayer(FFTLayer);
FFTCompressionLayer.setReferenceGrade(3);
FFTCompressionLayer.setCompressLevel(5);
% make the stack layer
functionStack = [sorReadLayer;sampleLayer;FFTLayer;FFTCompressionLayer];
FS = layerStack('SorghumFunctions',functionStack);
%% test node dependancy list
FS.getDependancyNodeList(FFTLayer.getUidType());
%% make network layers
% make the read layer
sorReadLayer1 = readLayer('sorghumStomataReader');
sorReadLayer1.setType('histogram');
% make the sample layer
sampleLayer1 = samplerLayer('SampleStomata');
% make and configure the fft layer
FFTLayer1 = fftLayer('sorghumStomataFFT');
FFTLayer1.setN(35);
% make and configure the compression layer
FFTCompressionLayer1 = compressionLayer(FFTLayer1);
FFTCompressionLayer1.setReferenceGrade(3);
FFTCompressionLayer1.setCompressLevel(5);
% make the stack layer
functionStack1 = [sorReadLayer1;sampleLayer1;FFTLayer1;FFTCompressionLayer1];
FS1 = layerStack('SorghumFunctions1',functionStack1);
%% make datapoint selector
dataS = dataPointSelector();
dataS.addAxis(sorFileList,'fileName');
dataS.addAxis(dataPointSelector.mat2cell(pointSet),'pointSet');
dataS.addAxis({opDomain},'domain');
dataS.addAxis({domainSize},'domainSize');
idx = dataS.rDrawN(10);
idx = dataS.rDrawN_fix1(10,10);
idx = dataS.rDrawN_fix1(2,1);
configureSample = dataS.drawToValue(idx);
FS.configure(configureSample);
FS1.configure(configureSample);
%% test single layer - read layer
clear res1
fluid = digitalFluid('test1',{inputPoint});
res1 = fluid.flow(sorReadLayer);
res1{1}
%% test single layer - layer stack
clear res1
fluid = digitalFluid('test1',{inputPoint});
res1 = fluid.flow(FS);
res1{1}
%% create the wide stack
wideStack = layerStack('finiteWideStack',[FS,FS1]);
%% persist test
fluid = digitalFluid('test1',{inputPoint});
topID = FFTLayer.bottomCE.uid;
fluid.persistData.toPersistFunction = @(x)strcmp(topID,x);
fluid.persistData.location = 'local';
res1 = fluid.flow(FS);
res1{1}
%% test the wide stack - with single input
clear res1
fluid = digitalFluid('test1',{inputPoint});
tic
res1 = fluid.flow(wideStack);
toc
res1{1}
%%
parpool(12)
%%
delete(gcp)
%% test the wide stack - with 10 inputs
clear res1
clear res2
fluid = digitalFluid('test1',inputPoints);
fluid.canParallel = false;
tic
res1 = fluid.flow(wideStack);
toc
size(res1)
res1{1}

fluid = digitalFluid('test1',inputPoints);
fluid.canParallel = true;
tic
res2 = fluid.flow(wideStack);
toc
size(res2)
res2{1}



fluid = digitalFluid('test1',TinputPoints(1:1000));
fluid.canParallel = false;
tic
res1 = fluid.flow(wideStack);
toc

fluid = digitalFluid('test1',TinputPoints(1:1000));
fluid.canParallel = true;
tic
res2 = fluid.flow(wideStack);
toc

%% test tall stack
clear res1
tic
fluid = digitalFluid('wholeTest',TinputPoints(1:10));
res1 = fluid.flow(FS);
size(res1)
res1{1}
toc
%% test with compute off and tracing only
clear res1
tic
fluid = digitalFluid('wholeTest',TinputPoints(1:1000));
fluid = digitalFluid('wholeTest',{inputPoint});
fluid.computeToggle(false);
res1 = fluid.flow(FS);
size(res1)
res1{1}
q = mksqlite('SELECT * FROM layerTable');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% where am i path number at node number make the graph

[nodeList,upid] = wideStack.topCE.getNodeList();
upid(1,:) = [];
connectionMatrix = zeros(numel(nodeList));

for e = 1:size(upid,1)
    childIDX = strcmp(nodeList,upid{e,1});
    parentIDX = strcmp(nodeList,upid{e,2});
    connectionMatrix(parentIDX,childIDX) = 1;
end
G = biograph(connectionMatrix,nodeList);
view(G);




%% reset database
mksqlite('open',':memory:');
mksqlite('CREATE TABLE layerTable (id INTEGER PRIMARY KEY AUTOINCREMENT,uid,pid,UNIQUE(uid,pid))');
mksqlite('CREATE TABLE currentLayer (uid)');
mksqlite(['INSERT INTO currentLayer (uid) VALUES (?)'],'null');
mksqlite('SELECT * FROM layerTable');

%% try build out over one image
% build out list
for e = 1:size(pointSet,1)
    inputList{e} = {sorFileList{1},pointSet(e,:),opDomain,domainSize};
end
inputList = inputList';
%% submit list
fluid = digitalFluid('test1',inputList(1:2));
fluid.flow(FS);
%%%%%%%%%%%%%%%%%%%%%%
%% test par for loop for handle objects
clear fluidTest
fluidTest = digitalFluid('test1',cell(100,1));
parfor e = 1:100
    %fluidTest.setData(e,e);
    fluidTest{e} = e;
end
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out X and Y for training - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
humanClickPath = '/mnt/tetra/nate/sorStomataClick/';
%humanClickPath = '/mnt/tetra/nate/maizeStomataClick/';
%humanClickPath = '/mnt/tetra/nate/seteriaStomataClick/';
trainFileList = {};
FileExt = {'mat'};
trainFileList = gdig(humanClickPath,trainFileList,FileExt,1);
dilateValue = 9;
%[labeledTrainingPackage] = generateMaskPackage(trainFileList,100,dilateValue,basisU,basisE);
[labeledTrainingPackage_CORN] = generateMaskPackage(trainFileList,100,dilateValue,basisU_CORN,basisE_CORN);

%%

flowData = traceNat('bug1');

flowData.flowDirection = 'f';
flowData.stopData.targetQueue = FFTLayer;

flowData.persistData.toPersist = false;
flowData.persistData.userName = 'nmiller';
flowData.persistData.location = 'irods';

flowData.jobData.uid = '1';


a = FS.flow(flowData,inputPoint{:});
%%
%% try running just the FFT layer
flowData = traceNat('bug1');

flowData.flowDirection = 'r';
flowData.stopData.targetQueue = FFTLayer;

flowData.persistData.toPersist = false;
flowData.persistData.userName = 'nmiller';
flowData.persistData.location = 'irods';

flowData.jobData.uid = '1';

[subIFFT,flowData] = FFTLayer.flow(flowData,inputPoint{:});

%%
inputPoint = {sorFileList{1},pointSet(1,:),opDomain,domainSize};

I = sorReadLayer.compute(inputPoint{:});

flowData.flowDirection = 'r';
flowData.history = [];
flowData.stopData.uid = sampleLayer.uid;
flowData.stopData.goFlag = false;
[subI,flowData] = sampleLayer.flow(flowData,inputPoint{:});

%% try running upto the compression layer
flowData.flowDirection = 'r';
flowData.history = [];

flowData.stopData.uid = FFTCompressionLayer.uid;
flowData.stopData.goFlag = false;

flowData.persistData.toPersist = false;
flowData.persistData.userName = 'nmiller';
flowData.persistData.location = 'irods';

flowData.jobData.uid = '1';


[subIFFT,flowData] = FFTCompressionLayer.flow(flowData,inputPoint{:});
%% try running the layer stack

flowData.flowDirection = 'r';
flowData.history = [];

flowData.stopData.uid = FS.uid;
flowData.stopData.goFlag = false;

flowData.persistData.toPersist = false;
flowData.persistData.userName = 'nmiller';
flowData.persistData.location = 'irods';

flowData.jobData.uid = '1';


[layerOutput,flowData] = FS.flow(flowData,inputPoint{:});




%%

dataS = dataPointSelector();
dataS.addAxis(sorFileList,'fileName');
dataS.addAxis(dataPointSelector.mat2cell(pointSet),'pointSet');
dataS.addAxis({opDomain},'domain');
dataS.addAxis({domainSize},'domainSize');

idx = dataS.rDrawN(10);
configurationSample = dataS.drawToValue(idx);

I = sorReadLayer.compute(configurationSample{1}{:});




%%



%% try running the extractionLayer
fileToRun = 1;
[pth,nm,ext] = fileparts(sorFileList{fileToRun});
customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
%customData{2} = globalIDX{sD};
AIextractionLayer(sorFileList{1},'stomata_histonormalize',{readerHistogram_sor},'stomata_fft',{35},domainData,subPointSet{sD},customData,'','');


%%

%TsorFileList = issueBulkTicket(sorFileList);
TmaizeFileList = issueBulkTicket(maizeFileList);
%TsetFileList = issueBulkTicket(setFileList);
rPath = '/iplant/home/nmiller/workIT/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_CORN2/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_SET/';
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_CORN_BETTER/';
[rPath,iticket] = issueTicket(rPath(1:end-1),10*numel(maizeFileList),'write');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = './outputFeatures/';
fibratedExtrationLayer = cFlow('extractionLayer');
fibratedExtrationLayer.setMCRversion('v930');
for e = 1:100
    for sD = 1%:numel(subPointSet)
        %[pth,nm,ext] = fileparts(TsorFileList{e});
        [pth,nm,ext] = fileparts(TmaizeFileList{e});
        %[pth,nm,ext] = fileparts(TsetFileList{e});
        customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
        customData{2} = globalIDX{sD};
        % for local run
        %fibratedExtrationLayer(TsorFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
        fibratedExtrationLayer(TmaizeFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
        %fibratedExtrationLayer(TsetFileList{e},'stomata_histonormalize',{toH},'stomata_fft',{35},domainData,subPointSet{sD},customData,oPath,rPath);
    end
    e
end
fibratedExtrationLayer.submitDag(auth,150,150);



%% JUNK BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create square grid - 
G = [];
g = [];
N = [50 60];        % this is 101 x 101 pixel square patch with 60 points along the edge
R = [20 40 60];

for r = 1:numel(R)
    [g(:,:,1),g(:,:,2)] = ndgrid(linspace(-R(r),R(r),N(1)),linspace(-R(r),R(r),N(1)));
    G(:,:,:,r) = cat(3,g(:,:,1),g(:,:,2));
end
%% create fraction grid to sample the images only at a percent (.35) of the image
% it was found to be better to sample many images at a smaller area
% than to sample few images over the whole area
imageFraction = .35;                % fraction of the image to sample
wholeSize = [512 512];              % whole size of the image
initSize = abs(min(G(:)));
[pt(:,:,1),pt(:,:,2)] = ndgrid(initSize:round(wholeSize(1)*imageFraction),initSize:round(wholeSize(2)*imageFraction));
ptSZ = size(pt);
pt = reshape(pt,[prod(ptSZ(1:2)) ptSZ(3)]);
%% generate tickets over both lists
TsorFileList = issueBulkTicket(sorFileList);
TmaizeFileList = issueBulkTicket(maizeFileList);
%% generate list for training
NI = 50;
trainFileList = {};
for e = 1:NI
    trainFileList{e} = sorFileList{e};
end
for e = 1:10
    trainFileList{end+1} = QClist{takeIDX(e)};
end
%% sample the fractional area for the training
for i = 1:numel(trainFileList)
    I = imread(trainFileList{i});
    sam = ba_interp2(I,pt(:,2),pt(:,1));
    sam = reshape(sam,ptSZ(1:2));
    imshow(sam,[]);
    title(num2str(i));
    drawnow
    samI(:,:,i) = sam;
end
%% collect clicks over the fractional area for training
for i = 51:numel(trainFileList)
    [col{i},row{i},~] = impixel(samI(:,:,i),[]);
end
%% make masks for the clicking data - to be used for training
MSK = [];
close all
for i = 1:NI
    msk = zeros(size(samI,1),size(samI,2));
    for p = 1:numel(col{i})
        msk(row{i}(p),col{i}(p)) = 1;
    end
    msk = imdilate(msk,strel('disk',5,0));
    out = flattenMaskOverlay(samI(:,:,i),logical(msk));
    imshow(out,[]);
    drawnow
    MSK = [MSK;msk(:)];
end