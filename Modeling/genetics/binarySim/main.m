%% test idea
BYTE_SIZE = 3;                              % log size of a byte 
MAX_B = 8;                                  % max bits to handle at once
%MAX_B = 12;                                 % max bits to handle at once
MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide
TALL = 10000;
xD = zeros(TALL,MAX_WID,'uint8'); % pre allocate data,'uint8'); % pre allocate data
%INT = uint8(randi((2^MAX_B)-1,TALL,1)); % int to test no zeros!! % generate random data
INT = double(randi((2^MAX_B)-1,TALL,1)); % int to test no zeros!! % generate random data
BYTE = floor(double(INT)*(2^BYTE_SIZE)^-1)+1;   
REM = rem(INT-1,2^BYTE_SIZE)+1;
POS = sub2ind(size(D),(1:size(D,1))',BYTE);
xD(POS) = bitset(D(POS),REM);

X = xD;
bit1 = randi(8,1);
bit2 = randi(8,1);
PROGRAM = randi(255,1,MAX_WID,'uint8');

%PROGRAM(1:end-1) = 0;

Y = uint8(bo(X,PROGRAM));

[STORE,DIS] = findProgram(10,2000,Y,X,10,MAX_WID,BYTE_SIZE,PROGRAM);


%%
mex -largeArrayDims -lpthread boo.cpp
%% inserted here - my sparse

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TALL = 100000;    
MAX_B = 8;
BYTE_SIZE = 3;                                                          % log size of a byte
MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);
TALL = 1000;


[INT] = samplePureStates(MAX_B,TALL);

%{
% for bi-linear logic
MAX_B = 4;
[INT] = samplePureStates(MAX_B,TALL);
[INT] = squeezeQbits([INT INT],4,3,'uint8');
MAX_B = 8;
%}
MAX_P_NUMBER = min(255,2^MAX_B);


br = bitRegister(MAX_B,INT);



%[X] = quantumEncode(MAX_B,INT,false);


[PROGRAM] = generateInitProgram(X,MAX_P_NUMBER,1,MAX_WID);

Y = evalProgram(PROGRAM,br);

evalProgramOutput = @(X)evalSimMetrics(X,Y,'hamming');

iterations = 200;
repeats = 10;
[STORE,DIS] = findProgram(1,[figure;figure],repeats,iterations,evalProgramOutput,br,Y,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER);
%% teach to add
Examples = 200;
X1 = randi(128,Examples,1,'uint8');
X2 = randi(128,Examples,1,'uint8');
tmpData = cat(2,X1,X2);
dataSum = [];
for b = 1:8
    dataSum(:,:,b) = squeezeQbits2(tmpData,{[b],[b]},3,'uint8');
end
br = bitRegister(2,dataSum);
PROGRAM = uint8([(2.^(0:7))*[[0 1 1 0 1 0 0 1]' [0 0 0 1 0 1 1 1]']]');

Y = evalProgram(PROGRAM,br);
Y = squeeze(Y(:,1,:));
Y = double(Y) * (2.^(0:7))';
X3 = X1 + X2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% 1) select program from 2^(2^b) as the grouping variable
% 2) generate N inputs
% 3) evaluate the program at the N inputs
% 4) 
%%%%%
%     
close all
%MAX_B = 2;                                  % the number of
MAX_B = 8;                                  % the number of
BYTE_SIZE = 3;                              % log size of a byte
MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);
TALL = 1000;
[INT] = samplePureStates(MAX_B,TALL);
[X] = quantumEncode(MAX_B,INT);
% generate program
MAX_P_NUMBER = min(255,2^MAX_B);
PROGRAM = randi(MAX_P_NUMBER,1,MAX_WID,'uint8');
PROGRAM(1) = 102;
PROGRAM(2:end) = 0;
% evaluate program with c-code
Y = boo(X,PROGRAM);
% generate phenotypes
DIS = [-1 .5 1 .5];
%DIS = [0 .5 0 .5];
UQ = unique(Y);
for g = 1:numel(DIS)/2
    IDX = (g-1)*2 + 1;
    gIDX = find(Y == UQ(g));
    PH(gIDX) = normrnd(DIS(IDX),DIS(IDX+1),numel(gIDX),1);
end
evalProgramOutput = @(X)evalSimMetrics(X,PH,'lda');
iterations = 200;
repeats = 10;



func = cFlow('findProgram');
func.setMCRversion('v930');
evalProgramOutput = @(X)evalSimMetrics(X,PH,'lda');
[a,b] = func(1,[],repeats,iterations,evalProgramOutput,X,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(10,10,auth);


[STORE,DIS] = findProgram(1,[figure;figure],repeats,iterations,evalProgramOutput,X,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM);
[J,sidx] = min(DIS(:,1));
classStructureBASE = boo(X,STORE(sidx,:));
[h,pval] = ttest2(PH(classStructureBASE==0),PH(classStructureBASE==1));
%%
func = cFlow('findProgram');
func.setMCRversion('v930');
evalProgramOutput = @(X)evalSimMetrics(X,PH,'lda');
[a,b] = func(1,[],repeats,iterations,PH,'lda',X,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(10,10,auth);



delta_h = [];
delta_pval = [];
classStructure = [];
for bit = 1:MAX_B
    value = bitget(INT,bit);
    tmpINT = bitset(INT,bit,~value);
    tmpX = quantumEncode(MAX_B,tmpINT);
    classStructure(:,bit) = boo(tmpX,STORE(sidx,:));
    [delta_h(bit),delta_pval(bit)] = ttest2(PH(classStructure(:,bit)==0),PH(classStructure(:,bit)==1));
    hd(bit) = bitwise_hamming(classStructureBASE,uint8(classStructure(:,bit)));
end
%%


%%
G = readtext('~/Downloads/genoCSV.csv');
P = readtext('~/Downloads/phenoCSV.csv');
M = readtext('~/Downloads/mapCSV.csv');
%%
CG = G(:,2:end);
CG = (strcmp(CG,'A'));
close all
MH = figure;
IH = figure;
disp = [figure;figure];
h = [];
pval = [];
MAX_B = 1;
MAX_P_NUMBER = min(255,2^MAX_B);


BYTE_SIZE = 3;
MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide


PH = cell2mat(P(2:end,4:end));
V = [];
CHECK = [];
for t = 120%2%1:size(CDCHECK,2)
    PH = phC(:,1);
    PH = cell2mat(P(2:end,4:end));
    PH = PH(:,t);
    PH = CDCHECK(:,t);
    h = [];
    pval = [];
    INTERACTION_VEC = zeros(size(CG,1),1);

    INTERACTION_VEC = zeros(size(CG,1));

    
    
    %for w = 1:(size(CG,1)^2)

    for w = 1:(size(CG,1)-MAX_B)

        
        
        %[wIDX(1) wIDX(2)] = ind2sub([size(CG,1) size(CG,1)],w);


       wIDX = w:w+(MAX_B-1);


        %wIDX = randi(size(CG,1),3,1);
        %wIDX = [201 3];


        INT = (2.^(0:(MAX_B-1))*CG(wIDX,:))';
        [X] = quantumEncode(MAX_B,INT);
        evalProgramOutput = @(X)evalSimMetrics(X,PH,'lda');
        iterations = 10;
        repeats = 5;
        br = bitRegister(MAX_B,INT);
%[STORE,DIS] = findProgram(1,[;],repeats,iterations,evalProgramOutput,br,Y,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER);

        [STORE,DIS] = findProgram(1,[],repeats,iterations,evalProgramOutput,br,PH,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER);
        [J,sidx] = min(DIS(:,1));
        classStructureBASE = boo(X,STORE(sidx,:));
        [h(w),pval(w)] = ttest2(PH(classStructureBASE==0),PH(classStructureBASE==1));


        STORE
        if  pval(w) < .001
            %{
            INTERACTION_VEC(wIDX(1),wIDX(2)) = pval(w);
            INTERACTION_VEC(wIDX(2),wIDX(1)) = pval(w);
            DIS = -log10(INTERACTION_VEC);
            DIS(size(CG,1)) = 0;
            DIS(isinf(DIS)) = 0;
            imshow(DIS,[]);
            drawnow
            w
            wIDX
%}
            %break
        end

        % simple
        [~,CHECK(w)] = ttest2(PH(INT==0),PH(INT==1));

        delta_h = [];
        delta_pval = [];
        classStructure = [];
        hd = [];
        for bit = 1:MAX_B
            value = bitget(INT,bit);
            tmpINT = bitset(INT,bit,~value);
            tmpX = quantumEncode(MAX_B,tmpINT);
            classStructure(:,bit) = boo(tmpX,STORE(sidx,:));
            [delta_h(bit),delta_pval(bit)] = ttest2(PH(classStructure(:,bit)==0),PH(classStructure(:,bit)==1));
            hd(bit) = bitwise_hamming(classStructureBASE,uint8(classStructure(:,bit)));
        end


        widx = find(hd > 10);
        if pval(w) < .01
            INTERACTION_VEC(wIDX(widx)) = INTERACTION_VEC(wIDX(widx)) + 1;
        end

        %{
        figure(MH);
        plot(-log10(pval))


        figure(IH);
        plot(INTERACTION_VEC)
        drawnow
%}
        t




        %figure;
        plot(-log10(CHECK))
        drawnow
    end

    V(:,t) = pval';



end
%%
%%
CMD = cell2mat(M(2:end,3));
fidx = [1;find(CMD==0);size(CMD,1)];
newC = [];
for e = 1:numel(fidx)-1
    tmp = cumsum(CMD(fidx(e):fidx(e+1)));
    tmp = tmp / max(tmp);
    tmp(1) = [];
    newC = [newC;tmp];
end

%%
%%
PH = cell2mat(P(2:end,4:end));
[phS phC phU phE phL phERR phLAM] = PCA_FIT_FULL(PH,3);
phG = cell2mat(P(2:end,2));
%%
% 59.153
PH = cell2mat(P(2:end,4:end));
UQ = unique(phG);
%PH = phC;
CDCHECK = [];
for tm = 1:size(PH,2)
    for u = 1:numel(UQ)
        %uidx = find(strcmp(phG,UQ{u}));
        uidx = phG == UQ(u);
        %uidx = phG
        tmpD = PH(uidx,tm);
        
        CDCHECK(u,tm) = mean(tmpD);
    end
end
%%
% simple
[~,CHECK(w)] = ttest2(PH(INT==0),PH(INT==1));

%%
[sweepD] = sweepPCA(phC,phE,phU,diag(phLAM).^.5,1:3,5);
%%
R = normrnd(0,20,1,5000);
EDGE = -128:128
[N,edges] = histcounts(R,EDGE);
Y = discretize(R,EDGE);
Y = uint8(Y);

close all
plot(N)
%%
cd  '/mnt/spaldingdata/nate/inMemCondor/compiledFunctions/mecka/DAG/737393692046938/functionInputs';
oPath = '/mnt/spaldingdata/nate/inMemCondor/compiledFunctions/mecka/DAG/smallDAG/functionInputs/';
%cd  '/mnt/spaldingdata/nate/inMemCondor/compiledFunctions/mecka/DAG/7373936945052654/functionInputs';
%oPath = '/mnt/spaldingdata/nate/inMemCondor/compiledFunctions/mecka/DAG/largeDAG/functionInputs/';

cdir = dir()
cdir(1:2) = [];
for e = 1:numel(cdir)
    fname = [cdir(e).folder filesep cdir(e).name];
    a = load(fname);
    [pth,nm,ext] = fileparts(a.arg2);
    newName = ['./' nm ext];
    fidx = strfind(newName,'#');
    newName = newName(1:(fidx(1)-1));
    a.arg2 = newName;
    a.arg5 = '';
    a.arg4 = './';
    [pth,nm,ext] = fileparts(fname);
    matName = [oPath nm ext];
    save(matName,'-struct','a')
    e
    numel(cdir)
end




