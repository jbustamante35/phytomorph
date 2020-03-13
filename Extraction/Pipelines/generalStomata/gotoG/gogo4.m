%% dig for training data
TrainFilePath = '/mnt/tetra/nate/RIL_train/';
TrainFileList = {};
FileExt = {'nms'};
TrainFileList = gdig(TrainFilePath,TrainFileList,FileExt,1);
%% gather click points
sampleSize = [250 250];
%c = {};
%r = {};
G = [];
D = [];
pSZ = [51 51];
hSZ = (pSZ-1)/2
[g1,g2] = ndgrid(-hSZ(1):hSZ(1),-hSZ(2):hSZ(2));
GR = [g1(:) g2(:)];
IMR = [];

%%%%
HOP = 4;
I = imread(TrainFileList{e});

    
    I = I(1:sampleSize(1),1:sampleSize(2));
g = im2colF(I,pSZ,[HOP HOP]);
%%%%

str = 1;
SKIP = size(g,2);
stp = str + SKIP - 1;
G = zeros([pSZ,size(g,2)*100]);
for e = 1:200%numel(TrainFileList)
    
    I = imread(TrainFileList{e});
    if e == 1
        J = I;
    end
    I = imhistmatch(I,J);
    
    
    
    
    I = I(1:sampleSize(1),1:sampleSize(2));
    %{
    if e > numel(c)
        [c{e},r{e},V] = impixel(I,[]);
    end
    %}
    
    
    
    g = im2colF(I,pSZ,[HOP HOP]);
    %g = bsxfun(@minus,g,min(g,[],1));
    %g = bsxfun(@mtimes,g,max(g,[],1));
    g = reshape(g,[pSZ size(g,2)]);
    %[dg1,dg2] = gradient(g);
    %g = (dg1.^2 + dg2.^2).^.5;
   
    %G = cat(3,G,g);
    G(:,:,str:stp) = g;
    str = stp + 1;
    stp = str + SKIP - 1;
    %{
    MASK = zeros(size(I));
    IDX = sub2ind(size(MASK),r{e},c{e});
    MASK(IDX) = 1;
    dMASK = imdilate(MASK,strel('disk',8));
    
    
    d = im2colF(dMASK,pSZ,[2 2]);
    %g = bsxfun(@minus,g,min(g,[],1));
    %g = bsxfun(@mtimes,g,max(g,[],1));
    d = reshape(d,[pSZ size(d,2)]);
    %[dg1,dg2] = gradient(g);
    %g = (dg1.^2 + dg2.^2).^.5;
    D = cat(3,D,d);
    %}
    
    IMR = [IMR;e*ones(size(d,3),1)];
    
    %out = flattenMaskOverlay(double(I),logical(dMASK));
    
    
    imshow(I,[]);
    drawnow
    e
end

%D = reshape(D,[size(D,1)*size(D,2) size(D,3)]);
G = reshape(G,[size(G,1)*size(G,2) size(G,3)]);
D = zeros(size(G));
%%
numPC = 200;
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL_T(G,numPC);
gS = reshape(gS,[51 51 size(gS,2)]);
%%
toMatch = 400;
DEEP = 3;
ZM = NaN*zeros(512,512);
MM = NaN*zeros(512,512);
[l1 l2] = ndgrid(linspace(hSZ(1)+1,512-hSZ(1)-1,11),linspace(hSZ(1)+1,512-hSZ(1)-1,11));
dl1 = round(diff(l1,1,1)*.5 + l1(1:end-1,:));
dl2 = round(diff(l2,1,2)*.5 + l2(:,1:end-1));
l1 = floor(l1);
l2 = floor(l2);
J1 = randi(10,size(l1,1),size(l2,2)) - 5;
J2 = randi(10,size(l1,1),size(l2,2)) - 5;
%l1 = J1 + l1;
%l2 = J2 + l2;
initI = randi(prod(size(l1)),1);
indexImageMask = zeros(size(l1));
indexImageMask(initI) = 1;
indexImage = zeros(size(l1));
randFunction = randi(size(G,2),1);
indexImage(initI) = randFunction;
toMatch = indexImage(find(indexImageMask));
IDX = sub2ind(size(ZM),GR(:,1)+l1(find(indexImageMask)),GR(:,2)+l2(find(indexImageMask)));

newFill = imdilate(indexImageMask,strel('diamond',1)) == 1 & indexImageMask == 0;
newDX = [];
oldDX = [];
newDX = [];
oldDX = [];
[newDX(:,1),newDX(:,2)] = find(newFill);
[oldDX(:,1),oldDX(:,2)] = find(indexImageMask);

ZM(IDX) = G(:,randFunction);
MM(IDX) = D(:,randFunction);
close all
imshow(ZM,[]);
%%
for loop = 1:20
    [ZM,MM,indexImageMask,indexImage,LABEL] = findMinTile(IMR,true,l1,l2,oldDX,toMatch,l1,l2,newDX,GR,G,D,DEEP,ZM,MM,indexImageMask,indexImage);
 
    newFill = imdilate(indexImageMask,strel('diamond',1)) == 1 & indexImageMask == 0;
    newDX = [];
    oldDX = [];
    [newDX(:,1),newDX(:,2)] = find(newFill);
    [oldDX(:,1),oldDX(:,2)] = find(indexImageMask);
    toMatch = indexImage(find(indexImageMask));
end

 ZMm = adapthisteq(ZM,'NumTiles',[15 15]);
 imshow(ZMm,[]);
%%
ZM2 = ZM;
indexImageMask2 = zeros(size(dl1));
indexImage2 = indexImageMask2;
newDX = [1 1];
DEEP = 1;
for loop = 1:20
    [ZM2,indexImageMask2,indexImage2] = findMinTile(false,l1,l2,oldDX,toMatch,dl1,dl2,newDX,GR,G,DEEP,ZM2,indexImageMask2,indexImage2);
 
    newFill = imdilate(indexImageMask2,strel('diamond',1)) == 1 & indexImageMask2 == 0;
    newDX = [];
    %oldDX = [];
    [newDX(:,1),newDX(:,2)] = find(newFill);
    %[oldDX(:,1),oldDX(:,2)] = find(indexImageMask2);
    %toMatch = indexImage(find(indexImageMask2));
end
 ZMm2 = adapthisteq(ZM2,'NumTiles',[15 15]);
 imshow(ZMm2,[]);
%%

%for r = 1:10
    for c = 1:10


        deltaShiftC = [0 45];
        totalShift = deltaShiftC*c;


        
        
        
        
        
        
        [~,sidx] = sort(delta);
        %imshow([reshape(G(:,sidx(DEEP)),[pSZ]) reshape(G(:,toMatch),[pSZ])],[]);

        Z1 = zeros(512,512);
        Z2 = zeros(512,512);
        IDX1 = sub2ind(size(Z1),GR(:,1)+hSZ(1)+1,GR(:,2)+hSZ(1)+1);
        IDX2 = sub2ind(size(Z1),GR(:,1)+hSZ(1)+1+totalShift(1),GR(:,2)+hSZ(1)+1+totalShift(2));
        Z1(IDX1) = G(:,toMatch);
        Z2(IDX2) = G(:,sidx(DEEP));
        Z1(Z1==0) = NaN;
        Z2(Z2==0) = NaN;
        if c == 1
            ZM = cat(3,ZM,Z1,Z2);
        else
            ZM = cat(3,ZM,Z2);
        end
        %

        ZM = nanmean(ZM,3);

        toMatch = sidx(DEEP);
        imshow(ZM,[]);
        drawnow

    end
%end
    imshow(ZM,[]);

%%
numPC = 15;


G = reshape(G,[size(G,1) size(G,2)*size(G,3)]);
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL_T(G,numPC);
gS = reshape(gS,[pSZ  size(G,2)/pSZ(2)]);
G = reshape(G,[pSZ size(G,2)/pSZ(2)]);
gC = reshape(gC,[numPC pSZ(1) size(gC,2)/pSZ(2)]);
OV(1) = 15;
BEGIN1 = gC(:,1:OV(1),:);
END1 = gC(:,end-OV(1)+1:end,:);
END1 = reshape(END1,[numPC*OV(1) size(END1,3)]);
BEGIN1 = reshape(BEGIN1,[numPC*OV(1) size(BEGIN1,3)]);



G = permute(G,[2 1 3]);
G = reshape(G,[size(G,1) size(G,2)*size(G,3)]);
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL_T(G,numPC);
gS = reshape(gS,[pSZ size(G,2)/pSZ(1)]);
G = reshape(G,[pSZ size(G,2)/pSZ(2)]);
gC = reshape(gC,[numPC pSZ(2) size(gC,2)/pSZ(1)]);
OV(2) = 15;
BEGIN2 = gC(:,1:OV(2),:);
END2 = gC(:,end-OV(2)+1:end,:);
END2 = reshape(END2,[numPC*OV(2) size(END2,3)]);
BEGIN2 = reshape(BEGIN2,[numPC*OV(2) size(BEGIN2,3)]);

G = permute(G,[2 1 3]);
%%
close all
toMatch = 9000;
MAX_ROW = 4;
MAX_COL = 10;
STRIP = G(:,:,toMatch);
stickEND_RIGHT_TOP = END1(:,toMatch);
stickEND_BOTTOM_TOP = END2(:,toMatch);
sitckEND_RIGHT_BOTTOM = NaN*zeros(size(stickEND_RIGHT_TOP));
FUNY = {};
for col = 1:MAX_COL
    FUNY{col} = sitckEND_RIGHT_BOTTOM;
end
Z = zeros(MAX_ROW*pSZ(1),MAX_COL*pSZ(2),1);
rowPTR = 1;
colPTR = 1;
Z(1:pSZ(1),1:pSZ(2)) = STRIP;
colNUM = 1;
rowNUM = 1;
DEEP = 10;
%OV = [1 1];
while rowNUM < MAX_ROW
    
    % match right
    if ~any(isnan(stickEND_RIGHT_TOP))
        delta = bsxfun(@minus,stickEND_RIGHT_TOP,BEGIN1);
        delta = sum(delta.*delta,1).^.5;
        
        if ~any(isnan(FUNY{colNUM}))
            here = 1;
            delta2 = bsxfun(@minus,FUNY{colNUM},BEGIN2);
            delta2 = sum(delta2.*delta2,1).^.5;
            delta = delta + delta2;
        end
        
        
        [~,sidx] = sort(delta);
        sel = sidx(DEEP);
        toMatch = sel;
        newTILE_RIGHT = G(:,:,sel);
        stickEND_RIGHT_TOP = END1(:,sel);
        tmp = END2(:,sel);
       
        
        Z(rowPTR:(rowPTR+pSZ(2)-1),(colPTR+pSZ(2)-1-OV(2)+1):(colPTR+pSZ(2)-1)) = ...
            .5*(newTILE_RIGHT(:,1:OV(1)) + Z(rowPTR:(rowPTR+pSZ(2)-1),(colPTR+pSZ(2)-1-OV(2)+1):(colPTR+pSZ(2)-1)));
        Z(rowPTR:(rowPTR+pSZ(2)-1),(colPTR+pSZ(2)):(colPTR+pSZ(2)+pSZ(2)-1-OV(1))) = newTILE_RIGHT(:,OV(1)+1:end);
    end
    
    % match down
    if ~any(isnan(stickEND_BOTTOM_TOP))
        
        delta = bsxfun(@minus,stickEND_BOTTOM_TOP,BEGIN2);
        delta = sum(delta.*delta,1).^.5;
        if ~any(isnan(sitckEND_RIGHT_BOTTOM))
            here = 1;
            delta2 = bsxfun(@minus,sitckEND_RIGHT_BOTTOM,BEGIN1);
            delta2 = sum(delta2.*delta2,1).^.5;
            delta = delta + delta2;
        end
        
        
        
        
        [~,sidx] = sort(delta);
        sel = sidx(DEEP);
        newTILE_DOWN = G(:,:,sel);
        
         
      
        
        Z((rowPTR+pSZ(2)-1-OV(2)+1):(rowPTR+pSZ(1)-1),colPTR:(colPTR+pSZ(1)-1)) = ...
            .5*(newTILE_DOWN(1:OV(2),:) + Z((rowPTR+pSZ(2)-1-OV(2)+1):(rowPTR+pSZ(1)-1),colPTR:(colPTR+pSZ(1)-1)));
        
        Z((rowPTR+pSZ(2)):(rowPTR+pSZ(2)+pSZ(2)-OV(2)-1),colPTR:(colPTR+pSZ(1)-1)) = newTILE_DOWN(OV(2)+1:end,:);
        
        
        if colNUM >= 2
            newFUNY{colNUM-1} = END2(:,sel);
        end
        
        sitckEND_RIGHT_BOTTOM = END1(:,sel);
        stickEND_BOTTOM_TOP = tmp;
        
        
        if colNUM == 1
            stickEND_BOTTOM_SAVE = END2(:,sel);
        end
        
        
        colNUM = colNUM + 1;
        if colNUM > (MAX_COL)
            
            
            % match down once more my good friend
             colPTR = colPTR + pSZ(2)-OV(1);
            delta = bsxfun(@minus,stickEND_BOTTOM_TOP,BEGIN2);
        delta = sum(delta.*delta,1).^.5;
        if ~any(isnan(sitckEND_RIGHT_BOTTOM))
            here = 1;
            delta2 = bsxfun(@minus,sitckEND_RIGHT_BOTTOM,BEGIN1);
            delta2 = sum(delta2.*delta2,1).^.5;
            delta = delta + delta2;
        end
        
        
        
        
        [~,sidx] = sort(delta);
        sel = sidx(DEEP);
        newTILE_DOWN = G(:,:,sel);
        
         
      
        
        Z((rowPTR+pSZ(2)-1-OV(2)+1):(rowPTR+pSZ(1)-1),colPTR:(colPTR+pSZ(1)-1)) = ...
            .5*(newTILE_DOWN(1:OV(2),:) + Z((rowPTR+pSZ(2)-1-OV(2)+1):(rowPTR+pSZ(1)-1),colPTR:(colPTR+pSZ(1)-1)));
        
        Z((rowPTR+pSZ(2)):(rowPTR+pSZ(2)+pSZ(2)-OV(2)-1),colPTR:(colPTR+pSZ(1)-1)) = newTILE_DOWN(OV(2)+1:end,:);
        
        if colNUM >= 2
            newFUNY{colNUM-1} = END2(:,sel);
        end
        %%%%%%
            
            colPTR = 1;
            rowPTR = rowPTR + pSZ(1) - OV(2);
            
            delta = bsxfun(@minus,stickEND_BOTTOM_SAVE,BEGIN2);
            delta = sum(delta.*delta,1).^.5;
            
            [~,sidx] = sort(delta);
            sel = sidx(1);
            newTILE_DOWN = G(:,:,sel);
        
            
            if (rowNUM+1) ~= MAX_ROW
            

                Z((rowPTR+pSZ(2)-1-OV(2)+1):(rowPTR+pSZ(1)-1),colPTR:(colPTR+pSZ(1)-1)) = ...
                    .5*(newTILE_DOWN(1:OV(2),:) + Z((rowPTR+pSZ(2)-1-OV(2)+1):(rowPTR+pSZ(1)-1),colPTR:(colPTR+pSZ(1)-1)));
 
                Z((rowPTR+pSZ(2)):(rowPTR+pSZ(2)+pSZ(2)-OV(2)-1),colPTR:(colPTR+pSZ(1)-1)) = newTILE_DOWN(OV(2)+1:end,:);
            end
            
            
            
            
            stickEND_RIGHT_TOP = END1(:,sel);
            stickEND_BOTTOM_TOP = END2(:,sel);
            sitckEND_RIGHT_BOTTOM = NaN*zeros(size(stickEND_RIGHT_TOP));
            imshow(Z,[]);
            waitforbuttonpress
            
            rowNUM = rowNUM + 1;
            colNUM = 1;
            rowPTR = rowPTR + pSZ(1) - OV(2);
            
            FUNY = newFUNY;
            
            
        else
            colPTR = colPTR + pSZ(2)-OV(1);
        end
        
        
       
        
        
    end
    
    
    
    imshow(Z,[]);
    waitforbuttonpress
    
    drawnow
end
%%
[regionDetectionDomain] = generateMultiResSampleDomain([800 800],[1],[20 20]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the data for each species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% return subset for John
subList = readtext('/home/nate/forJohn.csv');
subList{1}(1:3) = [];
nm = {};
for e = 1:numel(sorFileList)
    [~,nm{e}] = fileparts(sorFileList{e});  
    nm{e} = strrep(nm{e},'#','');
end
for e = 1:numel(subList)
    if strcmp(subList{e}(1:2),'16')
        subList{e}(1:2) = [];
    end
end
[~,~,IDX] = intersect(subList,nm);
oPath = '/mnt/tetra/nate/returnForJohn2/'
for e = 1:numel(IDX)
    I = imread(sorFileList{IDX(e)});
    filenameOUT = [oPath nm{IDX(e)} '.tif'];
    subList(e)
    imwrite(I,filenameOUT);
    
    
    imshow(I,[]);
    BOX = [40 40 size(I,2)-80 size(I,1)-80];
    rectangle('Position',BOX,'EdgeColor','r');
    saveas(gca,[oPath nm{IDX(e)} '_RED.tif']);
    
end
%%
e=3;
I = imread(sorFileList{IDX(e)});
imshow(I,[]);
BOX = [40 40 size(I,2)-80 size(I,1)-80];
rectangle('Position',BOX,'EdgeColor','r');

%% gather all maizeData
dataPath = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData%'];
dataPath = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData/2016RIL_good_quality%'];
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s) - issue write ticket(s) over the directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TmaizeFileList = issueBulkTicket(maizeFileList);
rPath = '/iplant/home/phytomorphuser/workITOUT_FAST_CORN_BETTER/';
[rPat, iticket] = issueTicket(rPath(1:end-1),10*numel(TmaizeFileList),'write');
%% gather all setariaData
dataPath = ['/iplant/home/leakey_cyverse/setariaData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        setFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
TsetFileList = issueBulkTicket(setFileList);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the data for each species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather the normalization data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toH = [];
for e = 1:20
    tmpH = imread(TmaizeFileList{e});
    toH = [toH;tmpH(:)];
    imshow(tmpH,[]);
    e
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scan for training data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = ['/iplant/home/phytomorphuser/workIT%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
fPatchFileList = {};
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.mat')
        fPatchFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
fPatchFileList
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull training set local
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%massDownload('/iplant/home/phytomorphuser/workIT/', '.mat','/mnt/tetra/nate/workIT/');
massDownload('/iplant/home/phytomorphuser/workITOUT_FAST_CORN2/', '.mat','/mnt/tetra/nate/workIT_CORN_TRAIN/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scan local for mat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMatFilePath = '/mnt/tetra/nate/workIT/';
cMatFilePath = '/mnt/tetra/nate/workIT_CORN_TRAIN/';
cMatFilePath = '/mnt/tetra/nate/workIT_SET_TRAIN/';
fPatchFileList = {};
FileExt = {'mat'};
fPatchFileList = gdig(cMatFilePath,fPatchFileList,FileExt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out basis vector - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[basisU,basisE] = buildGradeBasis(fPatchFileList,100);
%[basisU_CORN,basisE_CORN] = buildGradeBasis(fPatchFileList,100);
[basisU_SET,basisE_SET] = buildGradeBasis(fPatchFileList,100);
%[data] = myRemoteLoader(gradeFileList{1},'T');
%data = prepareData(fPatchFileList{1},basisU,basisE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% build out human click map(s) - from CONDOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a temp local location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpPath = ['/mnt/tetra/nate/tmp/' tempname filesep];
mkdir(tmpPath);
humanClickPath = '/mnt/tetra/nate/sorStomataClick/';
humanClickPath = '/mnt/tetra/nate/maizeStomataClick/';
humanClickPath = '/mnt/tetra/nate/maizeStomataClick/';
mkdir(humanClickPath);
toClick = 100;
for e = 1:toClick
    try
        tmpM = zeros(size(rawI));
        [pth,nm,ext] = fileparts(fPatchFileList{e});
        outFile = [humanClickPath nm '.mat'];
        fileName = fPatchFileList{e};
        if ~exist(outFile)
            fileList{1} = fileName;
            %fileList = xfer_get({fileName},tmpPath,0,0);
            % change from rI to rawI at UIUC - not sure why i renamed the
            % varaiable-  but i did
            load(fPatchFileList{e},'T','rawI','customData');

            tmpM(customData{2}) = 1;
            R = regionprops(logical(tmpM),'boundingBox');
            subI = imcrop(rawI,R(1).BoundingBox);
            uix = [];
            [uix(:,2),uix(:,1),~] = impixel(subI,[]);

            
            rI = subI;
            imshow(subI,[]);
            hold on
            plot(uix(:,2),uix(:,1),'r*')
            hold off
            drawnow
            save(outFile,'T','rI','customData','uix');
        else
            e
        end
    catch ME
        ME
    end
end
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
%% train labeled package - CONDOR workflow
%[AI_layer] = trainAIlayer(labeledTrainingPackage,[1 1]);
[AI_layer_CORN] = trainAIlayer(labeledTrainingPackage_CORN,[1 1]);
%% applyAI layer - CONDOR
[probMap] = applyAIlayer(fPatchFileList{1},AI_layer,basisU,basisE,[108 108]);
%% optimize AI output
opti_para = {};
%% optimize output layer with swarm and nelder mead
%[opti_para] = optimizeAIoutput(AI_layer,labeledTrainingPackage);
[opti_para_CORN] = optimizeAIoutput(AI_layer_CORN,labeledTrainingPackage_CORN,opti_para_CORN);
%% test opti
[probMap] = applyAIlayer(fPatchFileList{4},AI_layer,basisU,basisE,[108 108]);
[grade,ret1,ret2,BO] = opti(probMap,'',opti_para,labeledTrainingPackage.oI(1:108,1:108,4));
%%
applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,{'hello'},AI_layer,basisU,basisE,opti_para,'','');
%%
func = @(X)applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,{'hello'},AI_layer,basisU,basisE,opti_para,'./output/','');
applyFuncWrapper = partialFunction(func,'metaLayers_sorghum');
%applyFuncWrapper.publish();
%% load from web
webPath = './';
websave([webPath 'help.mat'],'https://de.cyverse.org/dl/d/A2E6EFA2-37C3-4EC6-9335-302C29F9B2CC/metaLayers_sorghum.mat');
load('./help.mat');
%% run over hand scored data
for e  = 135:numel(IDX)
    obj.func(sorFileList{IDX(e)});
end

%% compile for john
cdir = dir('/mnt/tetra/nate/output/');
cdir(1:2) = [];
OUT = {};
cnt = 1;
for e = 1:numel(cdir)
    [pth,nm,ext] = fileparts(cdir(e).name);
    if strcmp(ext,'.csv')
        toRead = ['/mnt/tetra/nate/output/' nm '.csv'];
        D = readtext(toRead);
        OUT{cnt,1} = nm(1:10);
        OUT{cnt,2} = size(D,1)-1;
        cnt = cnt +1;
    end
    e
end
cell2csv('/mnt/tetra/nate/output/compile.txt',OUT);
%%
GT = readtext('/mnt/tetra/nate/GT.csv');
%% model?
close all
D = readtext('/mnt/tetra/nate/output/EF0035_2_4_spotdata.csv');
%D = readtext('/mnt/tetra/nate/output/EF0054_1_1_spotdata.csv');
%D = readtext('/mnt/tetra/nate/output/EF0001_2_4_spotdata.csv');
csvFile = '/mnt/tetra/nate/output/EF0076_1_2_spotdata.csv';
D = readtext('/mnt/tetra/nate/output/EF0076_1_2_spotdata.csv');
I = imread(strrep(csvFile,'spotdata.csv','overlay.jpg'));
 h1 = figure;
        h2 = figure;
        OUT = {};
        cnt = 1;
        CNT = [];
for f = 1:numel(cdir)
    [pth,nm,ext] = fileparts(cdir(f).name);
    if strcmp(ext,'.csv')

        RNM{cnt} = nm;
        if ~contains(GT{cnt,1},nm(1:10))
            break
        end
        
        csvFile = ['/mnt/tetra/nate/output/' nm '.csv'];
        D = readtext(csvFile);
        I = imread(strrep(csvFile,'spotdata.csv','overlay.jpg'));
        hope = load(strrep(csvFile,'_spotdata.csv','.mat'),'pMAP');
        
        MSK = zeros(512-80);
        LOC = [];
        for e = 2:size(D,1)
            LOC(e-1,:) = cell2mat(D(e,2:3));
            LOC(e-1,:) = round(LOC(e-1,:));
            MSK(LOC(e-1,2),LOC(e-1,1)) = 1;
        end
        DIST = bwdist(MSK);
        uA = mean(cell2mat(D(2:end,1)));
        DIS = squareform(pdist(LOC));
        MM = [];
        for e = 1:size(LOC,1)
            tmp = DIS(e,:);
            tmp(e) = inf;
            [MM(e)] = min(tmp);
        end
        WHAT = mean(MM) + (uA).^.5
        WHAT = mean(MM)
        lambda = 1/mean(MM);
        PROB = exp(-lambda*DIST);
        %imshow(PROB,[]);
        better = ((512-80)^2/WHAT^2)
        orgin = size(D,1)-1
        mean([orgin better])


        rM = [];
       
        

        for e = 1:50
            rM(e) = mean(MM);
            [~,newLOC] = min(PROB(:));
            [newXY(2),newXY(1)] = ind2sub(size(PROB),newLOC);

            MSK = zeros(size(MSK));
            LOC = [LOC ; newXY];
            LOC = round(LOC);
            for i = 1:size(LOC,1)
                MSK(LOC(i,2),LOC(i,1)) = 1;
            end

            DIS = squareform(pdist(LOC));
            MM = [];
            for i = 1:size(LOC,1)
                tmp = DIS(i,:);
                tmp(i) = inf;
                [MM(i)] = min(tmp);
            end

            DIST = bwdist(MSK);
            lambda = 1/mean(MM);

            PROB = exp(-lambda*DIST);
%{
            figure(h1)
            imshow(I,[]);
            hold on
            plot(LOC(:,1),LOC(:,2),'g*')
            plot(LOC(end,1),LOC(end,2),'bo')
            drawnow
            hold on
            
            %}
                %waitforbuttonpress
            %imshow(PROB,[]);
            %drawnow
            %{
             figure(h2)
            plot(rM)
            drawnow
            waitforbuttonpress
            %}
            %waitforbuttonpress
            
            
           
        end
       
        
        
        [~,toC] = max(rM);
        GGG(cnt,:) = rM;
        
        plot(GGG(cnt,:));
        hold on
        plot(toC,GGG(cnt,toC),'ro');
        rec = GT{cnt,2} - (size(D,1)-1);
        plot(rec,GGG(cnt,rec),'go');
        waitforbuttonpress
        cnt = cnt + 1;
        hold off
        %{
        
        
        CNT(cnt,:) = [size(D,1)-1 size(D,1)-1 + toC];
        OUT{cnt,1} = nm;
        OUT{cnt,2} = size(D,1)-1;
        OUT{cnt,3} = size(D,1)-1 + toC;
        
        plot(cell2mat(GT(1:cnt,2)),CNT(:,2),'.');
        hold on
        plot(cell2mat(GT(1:cnt,2)),CNT(:,1),'r.');
        title([num2str(corr(cell2mat(GT(1:cnt,2)),CNT(:,1))) '--' num2str(corr(cell2mat(GT(1:cnt,2)),CNT(:,2)))])
        
        
        drawnow
        
        %}
        
        %{
        figure(h1);
        imshow(I,[]);
        hold on
        plot(LOC(1:size(D,1)-1 + toC,1),LOC(1:size(D,1)-1 + toC,2),'r*')
        drawnow
        %}
        %plot(CNT(:,1),CNT(:,2),'.')
        
        drawnow
       
        
    end
end

cell2csv('/mnt/tetra/nate/output/compile2.txt',OUT);



%%
func_CORN = @(X)applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU_CORN,basisE_CORN,AI_layer_CORN},domainData,pointSet,{'hello'},AI_layer_CORN,basisU_CORN,basisE_CORN,opti_para_CORN,'./output/','');
applyFuncWrapper = partialFunction(func_CORN,'metaLayers_corn');
applyFuncWrapper.publish();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate training data on CONDOR for sorghum data
% need to run through feature extraction first before clicking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate operational grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I = imread(sorFileList{1});
I = imread(maizeFileList{1});
%I = imread(setFileList{1});
imageSZ = size(I);
borderSize = [40 40];
[pointSet,pointSetSize] = genIXgrid(imageSZ,[1 1],borderSize);
innerSize = imageSZ - 2*borderSize;
numMini = [4 4];
miniSZ = ceil(innerSize.*numMini.^-1);
miniLabel = [floor((pointSet(:,1)-borderSize(1)-1).*miniSZ(1).^-1) floor((pointSet(:,2)-borderSize(2)-1).*miniSZ(2).^-1)];
[UQ,i1,i2] = unique(miniLabel,'rows');
subPointSet = {};
for u = 1:size(UQ,1)
    subPointSet{u} = pointSet(i2==u,:);
    globalIDX{u} = find(i2==u);
end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate sample domain-2 disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[X,Y] = ndgrid(-30:30,-30:30);
%opDomain{2} = [X(:) Y(:)];
%domainSize{2} = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% package input(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trimValues = R;
pointSetData = pointSet;
domainData{1} = opDomain;
domainData{2} = domainSize;
domainData{3} = trimValues;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% issue ticket(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run all layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = './outputFeatures/';
func = cFlow('cRunner');
func.setMCRversion('v930');
for e = 1:100%numel(TmaizeFileList)
    for sD = 1%:numel(subPointSet)
        % set file name
        [pth,nm,ext] = fileparts(maizeFileList{e});
        customData{1} = ['{name_' strrep(pth,filesep,'FILESEP') 'FILESEP' nm '}{patch_' num2str(sD) '}'];
        %customData{2} = globalIDX{sD};
        % for local run
        func(maizeFileList{e},oPath,rPath,customData);
       
        sD
    end
    e
end
func.submitDag(auth,400,400);
%%
%%%%%%%%%%%
% to publish generalized loader
applyFunc = @(X,oPath,rPath,OTHER)applyAllLayers(X,'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,OTHER,AI_layer,basisU,basisE,opti_para,oPath,rPath);
applyFuncWrapper = partialFunction(applyFunc,'sorghumStomataApply');
applyFuncWrapper.publish();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% single application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 1;
rtPath = '';
otPath = '';
[fdata] = applyAllLayers(sorFileList{e},'stomata_histonormalize',{toH},'fullApply',{35,basisU,basisE,AI_layer},domainData,pointSet,customData,AI_layer,basisU,basisE,opti_para,otPath,rtPath);
%%
cRunner(TsorFileList{e},oPath,rPath,customData);