%% this will find the json files for each v1 and v2
% v1 and v2 are different version of the algorithm on CyVerse
FilePath = '/mnt/tetra/nate/caliAnalysis/v1/';
FileList{1} = {};
FileExt = {'json'};
FileList{1} = fdig(FilePath,FileList{1},FileExt,0);
%
FilePath = '/mnt/tetra/nate/caliAnalysis/v2/';
FileList{2} = {};
FileExt = {'json'};
FileList{2} = fdig(FilePath,FileList{2},FileExt,1);
%% clip off the remote base and attach the loca
remoteBase = '/iplant/home/canibas/maizeData/seedlingData/';
localBase = '/home/nate/mntCyverse/seedlingData/';
%sudo ./mountCyverse canibas/maizeData/seedlingData/ /home/nate/mntCyverse/seedlingData/
%% try to gather the dates
basePath = '/home/nate/mntCyverse/seedlingData/';
cdir = dir(basePath);
cdir(1:2) = [];
for e = 1:numel(cdir)
    try
        cdir(e).num = datenum(cdir(e).name);
    catch ME
        cdir(e).num = NaN;
    end
end
S(1).name = 'B97 RILS';
S(1).startDate = '07-Jan-2018';
S(1).stopDate = '30-Jan-2018';

S(2).name = 'P39 RILS';
S(2).startDate = '06-Feb-2018';
S(2).stopDate = '28-Feb-2018';

S(3).name = 'WiDIV Rep1';
S(3).startDate = '15-Apr-2018';
S(3).stopDate = '03-May-2018';

S(4).name = 'WiDIV Rep2';
S(4).startDate = '09-May-2018';
S(4).stopDate = '30-May-2018';

for e = 1:numel(S)
    str = datenum(S(e).startDate);
    stp = datenum(S(e).stopDate);
    fidx = find(str <= [cdir.num] & [cdir.num] <= stp);
    tmp = cdir(fidx);
    tmpList = {};
    fileList = {};
    cnt = 1;
    for f = 1:numel(tmp)
        tmpList{f} = [basePath tmp(f).name filesep];
        tdir = dir(tmpList{f});
        tdir(1:2) = [];
        for l = 1:numel(tdir)
            fileList{cnt} = [tmpList{f} tdir(l).name];
            cnt = cnt + 1;
        end
    end
    
    S(e).folders = tmpList;
    S(e).fileList = fileList;
end
%%    
options = weboptions('Timeout',180);
webPath = '~';
websave([webPath 'maizeSeedlings_HOTFIX.mat'],'https://de.cyverse.org/dl/d/69389EB9-FDC1-4A64-BE4B-D975E2CE0634/maizeSeedlings_HOTFIX.mat',options);
load([webPath 'maizeSeedlings_HOTFIX.mat']);
funcToCall_ver1 = obj.func;

websave([webPath 'maizeSeedlings_V2.mat'],'https://de.cyverse.org/dl/d/0F3A630D-899E-4D6F-ACB3-CA46F253F1FD/maizeSeedlings_newCNN.mat',options);
load([webPath 'maizeSeedlings_V2.mat']);
funcToCall_ver2 = obj.func;   
%% test calls on a random image
testFile = S(1).fileList{1};
testFile = strrep(testFile,localBase,remoteBase);
funcToCall_ver1(testFile);
funcToCall_ver2(testFile);
%% issue ticket to test file
[testFileList] = issueBulkTicket({testFile});
%% pull out the inputs to the func1 and func2
testFile = testFileList{1};
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};

a = load([webPath 'maizeSeedlings_V2.mat']);
%a = load([webPath 'maizeSeedlings_HOTFIX.mat']);

arg1 = struct2cell(a.obj.editArgs);
func = func2str(a.obj.func);
p1 = strfind(func,'(');
p2 = strfind(func,')');
f1 = func((p2(1)+1):p1(2)-1);
cfunc1 = cFlow(f1);
%cfunc1.addDirectoryMap('output>/mnt/spaldingdata/nate/condorMapTest/');
cfunc1.setMCRversion('v930');
cfunc1(testFile,arg1{1},arg1{2},arg1{3},arg1{4},arg1{5},arg1{6},arg1{7},arg1{8},arg1{9});
%cfunc1(testFile,'./output/',[],arg1{1}{3},arg1{2},arg1{3},arg1{4},arg1{5},arg1{6});
cfunc1.submitDag(auth,50,50);
%% bulk issue over collections
for e = 1:numel(S)
    S(e).tFileList = {};
    for f = 1:numel(S(e).fileList)
        S(e).tFileList{f} = strrep(S(e).fileList{f},localBase,remoteBase);
    end
    S(e).tFileList = issueBulkTicket(S(e).tFileList);
end
%% pull out the inputs to the func1 and func2
collection = 3;
algoVersion = 2;
baseCollectionMapping = '/mnt/spaldingdata/nate/maizeSeedling/';
collectionMapping = [baseCollectionMapping 'collection' num2str(collection) filesep num2str(algoVersion) filesep];
mmkdir(collectionMapping);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if algoVersion == 1
    a = load([webPath 'maizeSeedlings_HOTFIX.mat']);
else
    a = load([webPath 'maizeSeedlings_V2.mat']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arg1 = struct2cell(a.obj.editArgs);
func = func2str(a.obj.func);
p1 = strfind(func,'(');
p2 = strfind(func,')');
f1 = func((p2(1)+1):p1(2)-1);
cfunc1 = cFlow(f1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfunc1.addDirectoryMap(['output>' collectionMapping]);
cfunc1.setMCRversion('v930');

for e = 1:numel(S(collection).tFileList)
    fileName = S(collection).tFileList{e};
    
    if algoVersion == 1
        cfunc1(fileName,'./output/',[],arg1{1}{3},arg1{2},arg1{3},arg1{4},arg1{5},arg1{6});
    else
        cfunc1(fileName,arg1{1},arg1{2},arg1{3},arg1{4},arg1{5},arg1{6},arg1{7},arg1{8},arg1{9},'./output/','',false);
    end
    
end
cfunc1.submitDag(auth,50,50);
%% dig for json files
collection = 4;
algoVersion = 2;
baseCollectionMapping = '/mnt/spaldingdata/nate/maizeSeedling/';
for c = 1:collection
    for a = 1:algoVersion
        dataPath = [baseCollectionMapping 'collection' num2str(c) filesep num2str(a) filesep];
        tic
        FilePath = dataPath;
        FileList = {};
        FileExt = {'json'};
        FileList = fdig(FilePath,FileList,FileExt,1);
        S(c).jsonList{a} = FileList;
        toc
    end
end
%% generate image file and json file db
%%%%%%%%%%%%%%%%%%%%%%%%%%
% open database
%mksqlite('open','/mnt/tetra/nate/caliAnalysis/dataBase.db');
mksqlite('open',':memory:');
%mksqlite('open','/mnt/scratch2/dataBase.db');
% set blob type
mksqlite('typedBLOBs', 1 );
% wrapping
mksqlite('param_wrapping', 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the fields based on json file name - select the first json file in
% the list
%jsonFile = FileList1{1};
jsonFile = S(1).jsonList{1}{3};
[pth,nm,ext] = fileparts(jsonFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% create image table
o1 = strfind(nm,'{');
o2 = strfind(nm,'}');
data = {'',imageFile,jsonFile};
FLDS = ['qrEQ,algoVersion,imageFile,jsonFile,'];
for f = 1:numel(o1)
    tmp = nm(o1(f)+1:o2(f)-1);
    sep = strfind(tmp,'_');
    key = tmp(1:sep-1);
    value = tmp(sep+1:end);
    data{end+1} = value;
    FLDS = [FLDS key ',']
end
FLDS(end) = [];
mksqlite(['CREATE TABLE image (id INTEGER PRIMARY KEY AUTOINCREMENT,' FLDS ')']);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% create plant table
% pull fields from json file
filetext = fileread(jsonFile);
jsonData = jsondecode(filetext);
flds = fields(jsonData.seedlingDoc(1));
FLDS = [];
for e = 1:numel(flds)
    FLDS = [FLDS flds{e} ','];
end
FLDS(end) = [];
insertS2 = repmat('?,',[1 numel(flds)+4]);
insertS2(end) = [];
insertS2 = ['(' insertS2 ')'];
mksqlite(['CREATE TABLE plant (id INTEGER PRIMARY KEY AUTOINCREMENT,image_id INTEGER,Genotype,imageFile,' FLDS ' ,FOREIGN KEY (image_id) REFERENCES image(id))']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for v = 1:numel(FileList)
    cnt = 1;
    for e = 1:numel(FileList{v})
        try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % json file name - read and decode
            jsonFile = FileList{v}{e};
            filetext = fileread(jsonFile);
            jsonData = jsondecode(filetext);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % image file name
            [pth,nm,ext] = fileparts(jsonFile);
            sidx = strfind(pth,filesep);
            imageFile = [pth(1:sidx(end)) nm(1:end-5) '.tiff'];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if the total count is less than the count_T
            % then get the QR box.
            if cnt < CNT_T
                I = imread(imageFile);
                [TqrBANK{e},qrCropBox(cnt,:)] = getQRcode_2(I,.35);
                cnt =  cnt + 1;
            else
                box = mean(qrCropBox,1);
                box(1:2) = box(1:2) - 100;
                box(3:4) = box(3:4) + 200;
                box = box*.35^-1;
                COL = [box(1) box(1)+box(3)];
                ROW = [box(2) box(2)+box(4)];
                I = imread(imageFile,'PixelRegion',{round(ROW) round(COL)});
                [TqrBANK{e},~] = getQRcode_2(I,.35);
            end
        catch ME
            ME
        end
    end
    qrBANK{v} = TqrBANK;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% need to check the QR codes on condor - too many images and remote otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% looks like I need to run this so that the code will gather the cropbox
% where the cnt is less than threshold
per = .35;
for collection = 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    baseCollectionMappingNames = '/mnt/spaldingdata/nate/maizeSeedlingNames_ver3/';
    mmkdir(baseCollectionMappingNames);
    collectionMapping = [baseCollectionMappingNames 'collection' num2str(collection) filesep];
    mmkdir(collectionMapping);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qrFunc = cFlow('getQRcode_2');
    qrFunc.setMCRversion('v930');
    qrFunc.addDirectoryMap(['output>' collectionMapping]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qrFunc.addSquidD('core-3.2.1.jar');
    qrFunc.addSquidD('javase-3.2.1.jar');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    iList = S(collection).tFileList;
    for e = 1:3%:numel(iList)
        %[localMSG{e}] = getQRcode_2(iList{e},per);
        qrFunc(iList{e},per);
    end
    
    qrFunc.submitDag(auth,50,50);
end
%% dig for name files 
FilePath = '/mnt/spaldingdata/nate/maizeSeedlingNames/';
nameFileList = {};
FileExt = {};
nameFileList = fdig(FilePath,nameFileList,FileExt,0);
%% test name checker
[~,testName] = fileparts(S(1).fileList{1});
checkName(testName,nameFileList);
getQRcode_2(S(1).tFileList{10},per);
%% populate db
cnt = 1;
checkON = true;
FAILS = {};
for v = 1:numel(FileList)
    for e = 1:numel(FileList{v})
        try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % json file name - read and decode
            jsonFile = FileList{v}{e};
            filetext = fileread(jsonFile);
            jsonData = jsondecode(filetext);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % image file name
            [pth,nm,ext] = fileparts(jsonFile);
            sidx = strfind(pth,filesep);
            imageFile = [pth(1:sidx(end)) nm(1:end-5) '.tiff'];
            oldQRdata = nm(1:end-5);
            nameBase = [pth filesep nm(1:end-5)];
            if checkON
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get QR data
                %I = imread(imageFile);
                %[newQRdata,qrCropBox] = getQRcode_2(I,.35);
                newQRdata = qrBANK{v}{e};
            else
                newQRdata = oldQRdata;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % is QR eq check
            if numel(newQRdata) == numel(oldQRdata)
                if all(oldQRdata == newQRdata)
                    isEQ = true;
                else
                    isEQ = false;
                end
            else
                isEQ = false;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % pull data from QR string
            toPull = newQRdata;
            o1 = strfind(toPull,'{');
            o2 = strfind(toPull,'}');
            data = {'',isEQ,v,imageFile,jsonFile};
            for f = 1:numel(o1)
                tmp = toPull(o1(f)+1:o2(f)-1);
                sep = strfind(tmp,'_');
                key = tmp(1:sep-1);
                value = tmp(sep+1:end);
                data{end+1} = value;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % insert image data
            insertS1 = repmat('?,',[1 numel(data)]);
            insertS1(end) = [];
            insertS1 = ['(' insertS1 ')'];
            m = mksqlite(['INSERT INTO image VALUES ' insertS1], data);



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop over plants
            genoConvert = {'A' 'B' 'C'};
            for p = 1:3
                croppedImageName = [nameBase strrep('{stage2croppedImage_%}.tif','%',num2str(p))];
                % select genotype from previous insert
                genoKey = ['Geno' genoConvert{p}];
                m = mksqlite(['SELECT ' genoKey ',Plot,Treatment FROM image WHERE id=' num2str(cnt)]);
                % set fixed data
                %data = {'',cnt,m.(genoKey),m.Plot,m.Treatment};
                data = {'',cnt,m.(genoKey),croppedImageName};
                for f = 1:numel(flds)
                    try
                        data{end+1} = jsonData.seedlingDoc(p).(flds{f}).data;
                        if iscell(data{end})
                            data{end} = NaN;
                        end
                    catch
                        data{end+1} = NaN;
                    end
                end
                mksqlite(['INSERT INTO plant VALUES ' insertS2], data);
            end
            
            
            
            e
            numel(FileList{v})

            cnt = cnt + 1;
        catch ME
            FAILS{end+1} = jsonFile;
            what = 1
        end
    end
end
%% number of datasets with fixed QR code
% check to make sure that the number of images makes sense
% do we want to load data from old QR - is there an joy here?
% started imaging at day 11 for cold - 11 on QR
% PCA explain to Cali
% dig into data "jumps" in raw data
% why are we missing some data - when we only have 1 curve
% A: second rep might have died
m = mksqlite(['SELECT qrEQ FROM image']);
badQR = sum([m.qrEQ] == 0)
%% get total genotypes
% got here on July 9, 2019
% tried running to make variables
% worked? -- answer = 
close all
cnt = 1;


clear Y
Y.bioMass = [];
Y.width = [];
Y.height = [];
Y.ID = [];

X = table;

 algoVersion = 1;
    dayBasis = 4:27;
  close all
treatmentTypes = {'Control','Cold'};
%treatmentTypes = {'Control'};
for treatment = 1:numel(treatmentTypes)
    treatement = treatmentTypes{treatment};
    
    genoList = mksqlite(['SELECT DISTINCT plant.Genotype FROM plant INNER JOIN image ON '...
        'image.id=plant.image_id '... 
        'WHERE image.Treatment=''' treatement  '''']);
  
   
    % for each genotype
    for toSelect = 1:numel(genoList)
        plotList = mksqlite(['SELECT DISTINCT image.Plot FROM image INNER JOIN plant ON '...
            'plant.image_id=image.id WHERE '...
            'image.Treatment=''' treatement  '''' ...
            'AND plant.Genotype=''' genoList(toSelect).Genotype '''' ...
            'AND image.algoVersion=' num2str(algoVersion)]);
        % for each plot
        for p = 1:numel(plotList)
            m = mksqlite(['SELECT image.imageFile,plant.id,image.PictureDay,plant.digitalBiomass,plant.plantHeight,plant.plantWidth FROM plant INNER JOIN image ON '...
                'image.id=plant.image_id WHERE' ... 
                 ' plant.Genotype=''' genoList(toSelect).Genotype  ... 
                 ''' AND image.Plot=''' plotList(p).Plot ...
                 ''' AND image.Treatment=''' treatement ...
                 ''' AND image.algoVersion=' num2str(algoVersion)]);
                
             
             
             
            bioMass = NaN*zeros(1,numel(dayBasis));
            height = NaN*zeros(1,numel(dayBasis));
            width = NaN*zeros(1,numel(dayBasis));
            ID = NaN*zeros(1,numel(dayBasis));
            % loop over days
            for d = 1:numel(m)
                m(d).PictureDay = str2num(m(d).PictureDay);
                
                % correct if is empty for is char
                if ischar(m(d).digitalBiomass) || isempty(m(d).digitalBiomass)
                    m(d).digitalBiomass = NaN;
                end
                
                % correct if is empty for is char
                if ischar(m(d).plantHeight) || isempty(m(d).plantHeight)
                    m(d).plantHeight = NaN;
                end
                
                % correct if is empty for is char
                if ischar(m(d).plantWidth) || isempty(m(d).plantWidth)
                    m(d).plantWidth = NaN;
                end
                
                fidx = find(dayBasis == m(d).PictureDay)
                bioMass(fidx) = m(d).digitalBiomass;
                height(fidx) = m(d).plantHeight;
                width(fidx) = m(d).plantWidth;
                ID(fidx) = m(d).id;
            end
            
            
            
             
            X(cnt,'Treatment') = {treatement};
            X(cnt,'Genotype') = {genoList(toSelect).Genotype};
            cnt = cnt + 1;
            Y.bioMass = [Y.bioMass;bioMass];
            Y.width = [Y.width;width];
            Y.height = [Y.height;height];
            Y.ID = [Y.ID;ID];
             
            
            
            plot(nanmean(Y.bioMass,1))
            %plot(Y.bioMass','r')
            drawnow
            
             %{
            % init vars for store
            bioMass = NaN*zeros(1,numel(dayBasis));
            ID = NaN*zeros(1,numel(dayBasis));

            % for each day
            for e = 1:numel(m)
                m1 = mksqlite(['SELECT id,PictureDay,imageFile,jsonFile,algoVersion FROM image WHERE id=' num2str(m(e).image) ...
                    ' AND algoVersion=' num2str(algoVersion) ...
                    ' AND Treatment=''' treatement '''' ...
                    ' AND Plot=''' plotList(p).Plot '''']);

                % if data is present
                if ~isempty(m1)
                    % get the index into the day vector
                    didx = find(dayBasis == str2num(m1.PictureDay));
                    % correct if is empty for is char
                    if ischar(m(e).digitalBiomass) || isempty(m(e).digitalBiomass)
                        m(e).digitalBiomass = NaN;
                    end
                    % assign to store
                    bioMass(1,didx) = m(e).digitalBiomass;
                    ID(1,didx) = m(e).id;
                end
            end
                 %}


             
             %{
            X(cnt,'Treatment') = {treatement};
            X(cnt,'Genotype') = {genoList(toSelect).Geno};
            cnt = cnt + 1;
            Y.bioMass = [Y.bioMass;bioMass];
            Y.ID = [Y.ID;ID];


            plot(nanmean(Y.bioMass,1))
            drawnow

%}
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BULK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fill in NaN
treatment = 'Cold';
fidx = find(strcmp(X.Treatment,treatment));
% fix biomass with nan
Y.bioMass = fixData(Y.bioMass);
Y.width = fixData(Y.width);
Y.height = fixData(Y.height);
%% PCA on data
tmpY = [Y.bioMass(:,1:14) Y.width(:,1:14) Y.height(:,1:14)];
kidx = find(~any(isnan(tmpY),2));
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(tmpY(kidx,:),2);
Y.bioMass_s = NaN*zeros(size(Y.bioMass));
Y.bioMass_s(kidx,1:14) = pS(:,1:14);

Y.width_s = NaN*zeros(size(Y.width));
Y.width_s(kidx,1:14) = pS(:,15:28);

Y.height_s = NaN*zeros(size(Y.height));
Y.height_s(kidx,1:14) = pS(:,29:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BULK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fill in NaN
basisIDX = [1:1+13];
treatment = 'Control';
fidx = find(strcmp(X.Treatment,treatment));
% fix biomass with nan
Y.bioMass(fidx,basisIDX) = fixData(Y.bioMass(fidx,basisIDX));
Y.width(fidx,basisIDX) = fixData(Y.width(fidx,basisIDX));
Y.height(fidx,basisIDX) = fixData(Y.height(fidx,basisIDX));
basisIDX = [8:8+13];
treatment = 'Cold';
fidx = find(strcmp(X.Treatment,treatment));
% fix biomass with nan
Y.bioMass(fidx,basisIDX) = fixData(Y.bioMass(fidx,basisIDX));
Y.width(fidx,basisIDX) = fixData(Y.width(fidx,basisIDX));
Y.height(fidx,basisIDX) = fixData(Y.height(fidx,basisIDX));
%% PCA on data

basisIDX_control = [1:1+12];
basisIDX_cold = [8:8+12];
treatment = 'Control';
fidxControl = find(strcmp(X.Treatment,treatment));
treatment = 'Cold';
fidxCold = find(strcmp(X.Treatment,treatment));
tmpY_control = [Y.bioMass(fidxControl,basisIDX_control) ...
    Y.width(fidxControl,basisIDX_control) ...
    Y.height(fidxControl,basisIDX_control)];
tmpY_cold = [Y.bioMass(fidxCold,basisIDX_cold) ...
    Y.width(fidxCold,basisIDX_cold) ...
    Y.height(fidxCold,basisIDX_cold)];
tmpY = [tmpY_control;tmpY_cold];
kidx = find(~any(isnan(tmpY),2));
% PCA decompose
[pS pC pU pE pL pERR pLAM] = PCA_FIT_FULL(tmpY(kidx,:),2);

bioMass_s = NaN*zeros(size(tmpY,1),numel(basisIDX_control));
width_s = NaN*zeros(size(tmpY,1),numel(basisIDX_control));
height_s = NaN*zeros(size(tmpY,1),numel(basisIDX_control));

bioMass_s(kidx,:) = pS(:,1:13);
width_s(kidx,:) = pS(:,14:26);
height_s(kidx,:) = pS(:,27:end);

Y.bioMass_s(fidxControl,basisIDX_control) = bioMass_s(1:numel(fidxControl),:);
Y.width_s(fidxControl,basisIDX_control) = width_s(1:numel(fidxControl),:);
Y.height_s(fidxControl,basisIDX_control) = height_s(1:numel(fidxControl),:);

Y.bioMass_s(fidxCold,basisIDX_cold) = bioMass_s(numel(fidxControl)+1:end,:);
Y.width_s(fidxCold,basisIDX_cold) = width_s(numel(fidxControl)+1:end,:);
Y.height_s(fidxCold,basisIDX_cold) = height_s(numel(fidxControl)+1:end,:);

Y.bioMass_s(fidxCold,1:basisIDX_cold(1)-1) = 0;
Y.width_s(fidxCold,1:basisIDX_cold(1)-1) = 0;
Y.height_s(fidxCold,1:basisIDX_cold(1)-1) = 0;

% set control and cold limits for fit
Y.bioMass_s(fidxControl,basisIDX_control(end)+1:end) = NaN;
Y.width_s(fidxControl,basisIDX_control(end)+1:end) = NaN;
Y.height_s(fidxControl,basisIDX_control(end)+1:end) = NaN;

% the below three lines are blocked out due to a change in the code on
% tina's bday
%{
Y.bioMass_s(fidxControl,basisIDX_cold(end)+1:end) = NaN;
Y.width_s(fidxControl,basisIDX_cold(end)+1:end) = NaN;
Y.height_s(fidxControl,basisIDX_cold(end)+1:end) = NaN;
%}
Y.bioMass_s(fidxCold,basisIDX_cold(end)+1:end) = NaN;
Y.width_s(fidxCold,basisIDX_cold(end)+1:end) = NaN;
Y.height_s(fidxCold,basisIDX_cold(end)+1:end) = NaN;

% set control and cold limits for raw
Y.bioMass(fidxControl,basisIDX_control(end)+1:end) = NaN;
Y.width(fidxControl,basisIDX_control(end)+1:end) = NaN;
Y.height(fidxControl,basisIDX_control(end)+1:end) = NaN;

% the below three lines are blocked out due to a change in the code on
% tina's bday
%{
Y.bioMass(fidxControl,basisIDX_cold(end)+1:end) = NaN;
Y.width(fidxControl,basisIDX_cold(end)+1:end) = NaN;
Y.height(fidxControl,basisIDX_cold(end)+1:end) = NaN;
%}
Y.bioMass(fidxCold,basisIDX_cold(end)+1:end) = NaN;
Y.width(fidxCold,basisIDX_cold(end)+1:end) = NaN;
Y.height(fidxCold,basisIDX_cold(end)+1:end) = NaN;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLOCK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sweep
[sweepD] = sweepPCA(pC,pE,pU,std(pC,1),1:2,10);
for p = 1:2
     close all
    h1 = figure;
    h2 = figure;
    h3 = figure;
    for s = 1:size(sweepD,2)
        t1 = squeeze(sweepD(p,s,1:13));
        t2 = squeeze(sweepD(p,s,14:26));
        t3 = squeeze(sweepD(p,s,27:end));
        figure(h1);
        plot(t1,'r');
        hold on
        figure(h2);
        plot(t2,'g');
        hold on
        figure(h3);
        plot(t3,'b');
        hold on
    end
    %waitforbuttonpress
end
%% scatter PCA
close all
plot(pC(numel(fidxControl)+1:end,1),pC(numel(fidxControl)+1:end,2),'b.')
hold on
plot(pC(1:numel(fidxControl),1),pC(1:numel(fidxControl),2),'r.')
fidx = find(all(pC(1:numel(fidxControl),:)<0,2));
fidx1 = find(pC(1:numel(fidxControl),1)>1*10^5 & pC(1:numel(fidxControl),2)>-.5*10^5);

sidx = 13;
plot(pC(fidx(sidx),1),pC(fidx(sidx),2),'go')
plot(pC(fidx1(sidx),1),pC(fidx1(sidx),2),'ko')
IDX = kidx(fidxControl(fidx(sidx)));
IDX1 = kidx(fidxControl(fidx1(sidx)));
figure;
plot(Y.bioMass(IDX,:),'g');
hold on
plot(Y.bioMass_s(IDX,:),'g--');
plot(Y.bioMass(IDX1,:),'k');
plot(Y.bioMass_s(IDX1,:),'k--');

%%
ids = Y.ID(IDX,:);
figure;
pause(1);
for m = 1:size(ids,2)
    try
        fileName = mksqlite(['SELECT imageFile FROM plant WHERE id=' num2str(ids(1,m))]);
        I = imread(fileName.imageFile);
        imshow(I,[]);
        drawnow
        %waitforbuttonpress
    catch
    end
end
%% plot unique genotypes
close all
UQ = unique(X.Genotype);
outSHEET = {};
cnt = 1;

for e = 1:numel(UQ)
    
   
    
    fidx = find(strcmp(X.Genotype,UQ{e}));
    subX = X(fidx,:);
    cUQ = unique(subX.Treatment);
    sig = Y.bioMass(fidx,:);
    ids = Y.ID(find(strcmp(X.Genotype,UQ{e})),:);
    sig_s = Y.bioMass_s(find(strcmp(X.Genotype,UQ{e})),:);
    
    for c = 1:numel(cUQ)
        gidx = find(strcmp(subX.Treatment,cUQ{c}));
        
        outSHEET{cnt,1} = UQ{e};
        outSHEET{cnt,2} = cUQ{c};
        U = nanmean(sig(gidx,:),1);
        for q = 1:numel(U)
            outSHEET{cnt,2+q} = U(q);
        end
        str = numel(U);
        U = nanmean(sig_s(gidx,:),1);
        for q = 1:numel(U)
            outSHEET{cnt,2+str+q} = U(q);
        end
        
        
        cnt = cnt + 1;
        
        if strcmp(cUQ{c},'Control')
            plot(nanmean(sig(gidx,:),1),'k');
            hold on
            plot(nanmean(sig_s(gidx,:),1),'k--');
            hold on
            plot(sig(gidx,:)','r')
            plot(sig_s(gidx,:)','m')
        else
            plot(nanmean(sig(gidx,:),1),'g');
            hold on
            plot(nanmean(sig_s(gidx,:),1),'g--');
            hold on
            plot(sig(gidx,:)','b')
            plot(sig_s(gidx,:)','c')
        end
    end
    
    
    
    hold off
    title(UQ{e});
    drawnow
    %waitforbuttonpress
    %{
    for m = 1:size(ids,2)
        try
            fileName = mksqlite(['SELECT imageFile FROM plant WHERE id=' num2str(ids(1,m))]);
            I = imread(fileName.imageFile);
            imshow(I,[]);
            drawnow
        catch
        end
    end
    %}
end
%%
oPath = '/mnt/tetra/nate/caliAnalysis/';
cell2csv([oPath 'cali_April_5_2019.csv'],outSHEET);
%%

%%
for e = 1:size(Y.bioMass,1)
    fidx = find(~isnan(Y.bioMass(e,:)));
    if any(fidx > 15)
        dayIdx = fidx(find(fidx > 15));
        ids = Y.ID(e,dayIdx);
        
        
        for i = 1:numel(ids)
            
            
            toCorrect = mksqlite(['SELECT image_id FROM plant WHERE id=' num2str(ids(i))]);
            toCorrect = mksqlite(['SELECT imageFile FROM image WHERE id=' num2str(toCorrect.image_id)]);
            fileToCorrect = toCorrect.imageFile;
            I = imread(fileToCorrect);
            [newQRdata,qrCropBox] = getQRcode_2(I,.35);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % pull data from QR string
            toPull = newQRdata;
            o1 = strfind(toPull,'{');
            o2 = strfind(toPull,'}');
            data = {'',isEQ,v,imageFile,jsonFile};
            for f = 1:numel(o1)
                tmp = toPull(o1(f)+1:o2(f)-1);
                sep = strfind(tmp,'_');
                key = tmp(1:sep-1);
                value = tmp(sep+1:end);
                data{end+1} = value;
            end
            
            
            
            
            
            data{end}
            dayBasis(dayIdx(i))
            fileToCorrect
            
            
            imshow(I,[]);
            waitforbuttonpress
        end
        
        
        
    end
end
%%
for e = 1:10
    
end
%%
m = mksqlite('SELECT * FROM plant')
%%

m = mksqlite('SELECT * FROM image')

            
            

