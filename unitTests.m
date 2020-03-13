%% unit tests for the re-structing of the phytoMorph pipeline
%% test storage blocks being glued to storage blocks
block1 = storageBlock({'hello'});
B = [block1;block1];
block2 = storageBlock(B);
block3 = [block1,block2];
%% repmatWide version
b1 = storageBlock({'hello'});
deepBlock = storageBlock([b1;copy(b1)]);
complexBlock = [copy(b1);deepBlock];
[newData] = storageBlock.repmatWide(complexBlock,3);
%% test if an array of inputBlocks are all inputBlocks
fprintf(['**************************************************************\n']);
fprintf(['Testing for homogeneous array types.\n']);
fprintf(['**************************************************************\n']);
in1 = inputBlock({'hello1'});
in2 = inputBlock({'hello2'});
inTotal = [in1,in2];
HomoTest = isa(inTotal,'inputBlock');
inHet = [block1,in1];
HetTest = isa(inHet,'inputBlock');
fprintf(['Test for homotype:' num2str(HomoTest) '\n']);
fprintf(['Test for heterotype:' num2str(HetTest) '\n']);
fprintf(['**************************************************************\n']);
%% recursive check for input data
inContainer = storageBlock([in1,in2]);
block1.containsInputData()
in1.containsInputData()
inContainer.containsInputData()
%% 
randNetwork = generateRandomNetwork();
expectedWholeOut = 0;
for e = 1:numel(randNetwork.layers)
    fprintf(['**************************************************************\n']);
    fprintf(['start TESTING LAYER NUMBER:' num2str(e) '\n']);
    fprintf(['**************************************************************\n']);
    singleLayerTest = randNetwork.layers(e);
    %singleLayerTest.view();
    fluid = digitalFluid('randTestFluid',{0});
    tic;fluid = fluid.flow(singleLayerTest);etim = toc;
    fprintf(['Flow test completed in:' num2str(etim) '(s).\n']);
    fprintf(['Output data is:' num2str(fluid.data{1}) '(s).\n']);
    fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
    fprintf(['Output data is sized:' num2str(size(fluid.data{1})) '\n']);
    newW = [];
    for e = 1:numel(expectedWholeOut)
        newW = [newW expectedWholeOut(e) + fluid.data{1}];
    end
    expectedWholeOut = newW;
    fprintf(['**************************************************************\n']);
    fprintf(['end TESTING LAYER NUMBER:' num2str(e) '\n']);
    fprintf(['**************************************************************\n']);
end
fprintf(['**************************************************************\n']);
fluid = digitalFluid('randTestFluid',{0});
tic;fluid = fluid.flow(randNetwork);etim = toc;
fprintf(['Flow test completed in:' num2str(etim) '(s).\n']);
fprintf(['Output data is:' num2str(fluid.data{1}) '(s).\n']);
fprintf(['Expected output is:' num2str(expectedWholeOut) '(s).\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['Output data is sized:' num2str(size(fluid.data{1})) '\n']);
%randNetwork.view();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make network layers
fprintf(['**************************************************************\n']);
fprintf(['Make network layers.Read->sample->fft->compression\n']);
fprintf(['**************************************************************\n']);
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
fprintf(['**************************************************************\n']);
%% make network layers
fprintf(['**************************************************************\n']);
fprintf(['Make network layers.Read->sample->fft->compression\n']);
fprintf(['**************************************************************\n']);
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
fprintf(['**************************************************************\n']);
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
fprintf(['**************************************************************\n']);
fprintf(['Test single read layer with single input.\n']);
fprintf(['Output should be cell array with first element as frozen tensor.\n']);
fprintf(['**************************************************************\n']);
fluid = digitalFluid('test1',inputPoint);
tic;fluid = fluid.flow(sorReadLayer);etim = toc;
fprintf(['Flow test completed in:' num2str(etim) '(s).\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['Output data is sized:' num2str(size(fluid.data{1})) '\n'])
fprintf(['**************************************************************\n']);
%% test single layer - layer stack
fprintf(['**************************************************************\n']);
fprintf(['Test layer stack upto and including compression layer for fft.\n']);
fprintf(['Output should be cell array with first element as compressed fft.\n']);
fprintf(['**************************************************************\n']);
clear res1
fluid = digitalFluid('test1',inputPoint);
tic;fluid = fluid.flow(FS);etim = toc;
fprintf(['Flow test completed in:' num2str(etim) '(s).\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['Output data is sized:' num2str(size(fluid.data{1})) '\n'])
fprintf(['**************************************************************\n']);
%% create the wide stack
wideStack = layerStack('finiteWideStack',[FS,FS1]);
%% test the wide stack - with single input
fprintf(['**************************************************************\n']);
fprintf(['Test wide layer with two branches..\n']);
fprintf(['**************************************************************\n']);
clear res1
fluid = digitalFluid('test1',inputPoint);
tic;res1 = fluid.flow(wideStack);etim = toc;
fprintf(['Flow test completed in:' num2str(etim) '(s).\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['Output data is sized:' num2str(size(res1.data{1})) '\n'])
fprintf(['**************************************************************\n']);
fprintf(['**************************************************************\n']);
%% persist test 1
fprintf(['**************************************************************\n']);
fprintf(['Test persist on read->sample->fft->compress.\n']);
fprintf(['Stop @ bottomCE fft.\n']);
fprintf(['**************************************************************\n']);
clear res1
fluid = digitalFluid('test1',inputPoint);
dData = persistQueue.generatePersistData('local','',1);
fluid.addPersistNode(FFTLayer.bottomCE.uid,dData);
fluid = fluid.flow(FS);
fprintf(['Fluid color is:' char(fluid.color) '.\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
tmp = load(fluid.data);
tmp = tmp.dFluid;
fprintf(['Loaded fluid color is:' char(tmp.color) '.\n']);
fprintf(['Output data is sized:' num2str(size(tmp.data{1})) '\n'])
% continue persist test 1 - continue flow
fprintf(['**************************************************************\n']);
fprintf(['Continue flow read->sample->fft->compress.\n']);
fprintf(['**************************************************************\n']);
fluid = fluid.flow(FS);
fprintf(['Fluid color is:' char(fluid.color) '.\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['Output data is sized:' num2str(size(fluid.data{1})) '\n'])
fprintf(['**************************************************************\n']);
fprintf(['**************************************************************\n']);
%% persist test 2 - branched persist
fprintf(['**************************************************************\n']);
fprintf(['Test persist on read->sample->fft->compress.\n']);
fprintf(['Stop @ bottomCE fft.\n']);
fprintf(['**************************************************************\n']);
fluid = digitalFluid('test1',inputPoint);
dData = persistQueue.generatePersistData('local','',1);
fluid.addPersistNode(FFTLayer.bottomCE.uid,dData);
fluid = fluid.flow(wideStack);
fprintf(['Fluid color is:' char(fluid.color) '.\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['**************************************************************\n']);
fprintf(['**************************************************************\n']);
tic;fluid = fluid.flow(wideStack);etim = toc;
fprintf(['Flow test completed in:' num2str(etim) '(s).\n']);
fprintf(['Manifold-q is sized:[' num2str(fluid.queueSize) '].\n']);
fprintf(['Output data is sized:' num2str(size(fluid.data{1})) '\n'])
fprintf(['**************************************************************\n']);
fprintf(['**************************************************************\n']);



%%
unit = @(x)x+1;
unit(1)
N = randi(5);

wholeStack = [];

for w = 1:N

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate a layerStack
    H = randi(3);
    tmpStack = [];
    for h = 1:H
        layerName = ['randLayer_' num2str(w) '_' num2str(h)];
        layer = functionHandleLayer(layerName,unit);
        tmpStack = [tmpStack;layer];
    end
    layerName = ['branchStack_' num2str(w)];
    tmpStack = layerStack(layerName,tmpStack);



    wholeStack = [wholeStack , tmpStack];
end

layerName = ['Total Stack'];
totalStack = layerStack(layerName,wholeStack);
totalStack.view()
%%

%%


    blockA = storageBlock({'A'});
    blockB = storageBlock(repmat(blockA,[1 3]));
    extractBlock = blockB.data(2);
    test = storageBlock(repmat(extractBlock,[1 2]));
    test.data(1).data{1} = 'B';

    inputBl


    blockA = storageBlock({'A'});
    blockB = storageBlock(repmat(blockA,[1 3]));
    extractBlock = copy(blockB.data(2));
    test = storageBlock(repmat(extractBlock,[1 2]));
    test.data(1).data{1} = 'B';
    
    clear A
    A = storageBlock({'A'});
    A = myModTest(A,4);

%}