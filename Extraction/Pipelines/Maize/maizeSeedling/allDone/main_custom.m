dataPath = ['/iplant/home/hirsc213/maizeData/seedlingData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nef')
        wholeFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
[wholeFileListT] = issueBulkTicket(wholeFileList);
%% cali scan
dataPath = ['/iplant/home/canibas/maizeData/seedlingData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nef')
        cali_wholeFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%% issue tickets over cali
[Tcali_wholeFileList] = issueBulkTicket(cali_wholeFileList);
%% create drop zone for data to land
rPath = '/iplant/home/phytomorphuser/cali_March_21_2019/';
system(['imkdir ' rPath]);
[rPath,iticket] = issueTicket(rPath(1:end-1),10*numel(Tcali_wholeFileList),'write');
%% try to re-run cali list with the partial function version
func = @(F,T)smartMain_new_ver1(F,coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,GMModel,GMModel2,E,U,GMModel_Color,'./output/',T,false);
pfa = partialFunction(func,'Cali_Seedling_March_21_2019');
pfa.publish();

ggC = cFlow('myPFwrapper');
ggC.setMCRversion('v930');
ggC.setMemory('8000');
res = {};
for i = 1:2%numel(TsorFileList)
   res{i} = ggC('Cali_Seedling_March_21_2019',Tcali_wholeFileList{i},rPath);
   fprintf(['Done building job:' num2str(i) '\n'])
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
ggC.submitDag(auth,150,150);

%%
oPath = './output/';
numJobs = numel(wholeFileListT);
[remoteOutputLocation,iticket] = issueTicket('/iplant/home/phytomorphuser/maizeSeedlingsLargeReturn',5*numJobs,'write');
rPath = remoteOutputLocation;
%%
func = cFlow('smartMain_new_ver1');
func.setMCRversion('v930');
func.setMemory('8000');

for e = 1:1000%numel(wholeFileListT)
    func(wholeFileListT{e},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,GMModel,GMModel2,E,U,GMModel_Color,oPath,rPath,false);
end
func.submitDag(auth,150,150);
%% publish again
%%%%%%%%%%%
% to publish smartMain
% run this with variables created from main
func = @(X)smartMain_new_ver1(X,coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,GMModel,GMModel2,E,U,GMModel_Color,'./output/','',false);
pF = partialFunction(func,'maizeSeedlings_newCNN');
pF.publish();

%%
smartMain_new_ver1(wholeFileListT{93},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,GMModel,GMModel2,E,U,GMModel_Color,oPath,rPath,false);
%%
load('~/junkList.mat','subFileList');
%%
cropData = getSeedlingCropBoxes(subFileList{1},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine);
%%
func = cFlow('getSeedlingCropBoxes');
func.setMCRversion('v930');
for e = 1:numel(subFileList)
    cropBoxData{e} = func(subFileList{e},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine);
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,150,150);
%%
for e = 1:numel(cropBoxData)
    d{e} = cFlowLoader(cropBoxData{e});
end
%% look through results from condor
FilePath = '/mnt/tetra/nate/maizeSeedlings2/';
FileList_cr = {};
FileExt = {'tif'};
FileList_cr = gdig(FilePath,FileList_cr,FileExt,1);
sigStack = [];
disp = true;
CLS = [];
for e = 1:numel(FileList_cr)
    if contains(FileList_cr{e},'croppedImage')
        Io = imread(FileList_cr{e});
        I = rgb2hsv(Io);
        I = I(:,:,3);
        sig = mean(I,2);
        sig = sig((end-200):end);
        sig = gradient(sig);
        
        sig = im2col(sig,[7 1],'sliding');
        if disp
            tC = PCA_REPROJ_T(sig,E,U);
            tC = tC(1:180);
            
            [~,midx] = min(tC);
            
            midx = size(I,1) - (200 - midx);
            
            
            Ic = Io(1:midx,:,:);
            Lab = rgb2lab(Ic);
            sz = size(Lab);
            CL = reshape(Lab(:,:,2:3),[prod(sz(1:2)) 2]);
            
            
            
            kidx = GMModel_Color.cluster(CL);
            kidx = reshape(kidx,sz(1:2));
            kidx = bwlarge(kidx==1);
            
            
            
            CL = CL(randperm(size(CL,1)),:);
            
            CLS = [CLS;CL(1:5000,:)];
            
            
            imshow(I,[]);
            hold on
            plot(1:size(I,2),midx*ones(1,size(I,2)),'r');
            drawnow
            hold off
        end
        sigStack = [sigStack sig];
        %plot(sig)
        drawnow
    end
end
%% test run mask 
for e = 1:numel(FileList_cr)
    if contains(FileList_cr{e},'croppedImage')
        Io = imread(FileList_cr{e});
        [kidx] = getMaskFromCropped(Io,E,U,GMModel_Color);
        imshow(kidx,[]);
        drawnow
    end
end
%%
GMModel_Color = fitgmdist(CLS,2);
%%
[S C U E L ERR LAM] = PCA_FIT_FULL_T(sigStack,1);
%%
SAM = [];
SAM2 = [];
for e = 1:numel(d)
    for p = 1:numel(d{e}.croppedImage)
        img = d{e}.croppedImage{p}/255;
        %img = imfilter(img,fspecial('gaussian',[31 31],5),'replicate');
        
        LAB = rgb2lab(img);
        sz = size(LAB);
        sam = LAB(:,:,1:3);
        sam = reshape(sam,[prod(sz(1:2)) 3]);
        
        if exist('GMModel')
            kidx = GMModel.cluster(sam);
            kidx = reshape(kidx,sz(1:2));
            plantMask = kidx == 1;
            plantMask = bwareaopen(plantMask,500);
            fidx = find(plantMask);
            %out = flattenMaskOverlay(d{e}.croppedImage{p}/255,plantMask);
          
            
            kidx2 = GMModel2.cluster(sam(fidx,1:2));
            kidx(fidx) = kidx2+3;
            out = kidx==4;
            
            plantMask = out;
            plantMask = bwareaopen(plantMask,400);
            plantMask = imopen(plantMask,strel('disk',3,0));
            plantMask = imclearborder(plantMask);
            plantMask = bwareaopen(plantMask,400);
            [plantMask] = connectPlant2(plantMask);
            
            
            
            plantMaskSub = imcrop(plantMask);
           
            %plantMaskSub = logical(interp2(plantMaskSub,2));
            
            outerMask = (bwdist(plantMaskSub));
            outerMask = (outerMask > 0).*(outerMask - 1);
            innerMask = (bwdist(~plantMaskSub));
            distT = innerMask - outerMask;
            E = edge(plantMaskSub);
            SK = bwmorph(plantMaskSub,'skeleton',inf);
            dE = bwdist(E);
            distT = dE.*plantMaskSub - dE.*(~plantMaskSub);
            %distT = imfilter(distT,fspecial('gaussian',[31 31],5),'replicate');
            th = .1;
            mlt = .2;
            funcS = @(X)(1 * (1 + exp(-.5*(X - 7))).^-1);
            
            for r = 1:5000
                MLT = del2(distT);
                newD = mlt*funcS(distT).*MLT;
                %newD = imfilter(newD,fspecial('gaussian',[31 31],5),'replicate');
                %newD = mlt*MLT;
               
                %distT = imfilter(distT,fspecial('gaussian',[31 31],5),'replicate');
                dM = distT > th*max(distT(:));
                %imshow(-del2(distT) > .1,[]);
                %imshow(distT.*dM,[]);
                %F1 = distT < (1+th) & distT > (1-th);
                %F2 = MLT < .001;
                %F3 = newD < -.005 & F1;
                
                
                plot(distT(find(plantMaskSub)),newD(find(plantMaskSub)),'.');
                hold on
                plot(distT(find(SK)),newD(find(SK)),'r.');
                hold off
                distT = distT + newD;
                %{
                R = regionprops(F3,F3.*newD,'MaxIntensity',');
                for t = 1:numel(R)
                    
                end
                %}
                %imshow(distT.*dM,[])
                %contour(distT,linspace(0,10,10));
                %ksdensity(newD(:))
                title([num2str(r) '--' num2str(max(distT(:)))])
                drawnow
            end
            
            plantSkeleton = bwmorph(plantMask,'skeleton',inf);
            endPoints = bwmorph(plantSkeleton,'endpoints');
            branchPoints = bwmorph(plantSkeleton,'branchpoints');
            skelpoint = [];
            [skelpoint(:,2),skelpoint(:,1)] = find(plantSkeleton);
            endpoint = [];
            [endpoint(:,2),endpoint(:,1)] = find(endPoints);
            branchpoint = [];
            [branchpoint(:,2),branchpoint(:,1)] = find(branchPoints);
            
            A = Radjacency(skelpoint', 2^.5+eps);
            EP = [];
            for r = 1:size(endpoint,1)
                EP(r) = find(all(bsxfun(@eq,skelpoint,endpoint(r,:)),2));
            end
            BP = [];
            for r = 1:size(branchpoint,1)
                BP(r) = find(all(bsxfun(@eq,skelpoint,branchpoint(r,:)),2));
            end
            
            G = graph(A);
            
            K = shortestpathtree(G,EP(1),EP(2));
            
            
            
            
            out = flattenMaskOverlay(img,logical(plantMask));
            
            imshow(out,[]);
            hold on
            plot(endpoint(:,1),endpoint(:,2),'r.');
            plot(branchpoint(:,1),branchpoint(:,2),'g.');
            plot(endpoint(1,1),endpoint(1,2),'r*');
            plot(endpoint(2,1),endpoint(2,2),'r*');
            plot(skelpoint([1 2 5 9],1),skelpoint([1 2 5 9],2),'b');
            drawnow
            
            waitforbuttonpress
            
            
        else
            out = LAB(:,:,3);
            fidx = fidx(randperm(numel(fidx)));
            SAM2 = [SAM2;sam(fidx(1:2000),:)];
            sam = sam(randperm(size(sam,1)),:);
            SAM = [SAM;sam(1:5000,:)];
            imshow(out,[]);
            drawnow
            waitforbuttonpress
        end
        
       
    end
    e
end
%%
GMModel = fitgmdist(SAM,2);
%%
GMModel2 = fitgmdist(SAM2(:,1:2),3);
%% run local on smart main for modeling
smartMain_new_ver1(FileList{1},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,GMModel,GMModel2,E,U,GMModel_Color,oPath,rPath,false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% local focus
FilePath = '/mnt/tetra/nate/seedlingDATApile/maizeData/seedlingData/';
FilePath = '/mnt/tetra/nate/hirschSampleRAWWHOLE/';
FilePath = '/mnt/tetra/nate/caliSampleRAWWHOLE/';
FilePath = '/mnt/tetra/nate/projectData/maizeSeedling/NEFs/';
FileList = {};
FileExt = {'nef'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert 5000 images to tiff
trainPath = '/mnt/tetra/nate/projectData/maizeSeedling/NEF_to_tiff/MN/';
FileList = FileList(randperm(numel(FileList)));
parfor e = 1:numel(FileList)
    fprintf(['Starting image conversion\n']);tic
    [pth,nm,ext] = fileparts(FileList{e});
    I = imread(FileList{e});
    imwrite(I,[trainPath nm '.tif']);
    fprintf(['Ending image conversion in:' num2str(toc) '\n']);
end
%% convert 5000 images to tiff
trainPath = '/mnt/tetra/nate/projectData/maizeSeedling/NEF_to_tiff/UW/';
FileList = FileList(randperm(numel(FileList)));
parfor e = 1:numel(FileList)
    try
        fprintf(['Starting image conversion\n']);tic
        [pth,nm,ext] = fileparts(FileList{e});
        newName = [trainPath nm '.tif'];
        if ~exist(newName,'file')
            I = imread(FileList{e});
            imwrite(I,newName);
        end
        fprintf(['Ending image conversion in:' num2str(toc) '\n']);
    catch ME
        getReport(ME)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD NETWORK FOR RED CORNER DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gather training data for the red corners
% sub task crop QR and build crop table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
FilePath = '/mnt/tetra/nate/projectData/maizeSeedling/NEF_to_tiff/';
FileList = {};
FileExt = {'tif'};
FileList = gdig(FilePath,FileList,FileExt,1);
qrPath = '/mnt/tetra/nate/projectData/maizeSeedling/QR_mixed/';
QRTable = table;
FileList = FileList(randperm(numel(FileList)));
deltaTM = [];

parfor e = 1:numel(FileList)
    redPOINTS = [];
    try
    
        fprintf(['Starting QR data gather for:' num2str(e) ':' num2str(numel(FileList)) '\n']);tic
        [pth,nm,ext] = fileparts(FileList{e});
        qrFileName = [qrPath nm ext];

        if ~exist(qrFileName,'file')

            I = imread(FileList{e});


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% start QR data gather
            Lab = rgb2lab(I);
            RED = Lab(:,:,2) > 25.5;
            RED = imclose(RED,strel('square',51));
            RED = bwlarge(RED);
            R = regionprops(RED);
            box = R(1).BoundingBox;
            box(1:2) = box(1:2) - pad;
            box(3:4) = box(3:4) + 2*pad;
            frame = imcrop(RED,box);
            frame = imclose(frame,strel('square',11));
            holes = imfill(frame,'holes') - frame;
            smallHoles = holes - bwareaopen(holes,1000);
            frame = logical(frame + smallHoles);
            qr = imcrop(I,box);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Lab = rgb2lab(qr);
            BLUE = Lab(:,:,3) < 0;
            BLUE = bwareaopen(BLUE,500);
            BLUE = imclose(BLUE,strel('square',13));
            BLUE = imclearborder(BLUE);
            BLUE = bwlarge(BLUE,24);
            %BLUEframe = imcrop(BLUE,box);
            BLUEframe = BLUE;
            BLUEframe = imdilate(BLUEframe,strel('square',9));
            %BLUEframe = imdilate(BLUEframe,strel('line',13,0));
            %BLUEframe = imdilate(BLUEframe,strel('line',13,90));
            %BLUEframe = imerode(BLUEframe,strel('line',13,0));
            %BLUEframe = imerode(BLUEframe,strel('line',13,90));
            %BLUEframe = imerode(BLUEframe,strel('square',5));
            holes = imfill(BLUEframe,'holes') - BLUEframe;
            smallHoles = holes - bwareaopen(holes,1000);
            BLUEframe = logical(BLUEframe + smallHoles);
            BLUEskeleton = bwmorph(BLUEframe,'skel',inf');
            BLUEskeleton = bwmorph(BLUEskeleton,'spur',inf);
            BLUEskeleton = bwmorph(BLUEskeleton,'skel',inf');
            BLUEskeleton = imdilate(BLUEskeleton,strel('square',1));
            BLUEsquares = imfill(BLUEskeleton,'holes');
            BLUEsquares = imdilate(BLUEsquares,strel('square',4));
            dB = bwboundaries(BLUEsquares);
            boxCorners = [];
            for b = 1:numel(dB)
                K = cwtK_closed_peaks(dB{b}(1:(end-1),:),{11},50,1.2*10^-8);
                boxCorners(:,:,b) = flipdim(dB{b}(K.pIDX,:),2);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}




            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            skeleton = bwmorph(frame,'skel',inf');
            skeleton = bwmorph(skeleton,'spur',inf);
            skeleton = bwmorph(skeleton,'skel',inf');
            branchPoints = bwmorph(skeleton,'branchpoints',1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% find the point orders
            br = [];
            [br(:,2),br(:,1)] = find(branchPoints==1);
            [~,TOPPOINT_IDX] = min(br(:,2));
            [~,LEFTPOINT_IDX] = min(br(:,1));
            [~,RIGHTPOINT_IDX] = max(br(:,1));
            MIDDLEPOINT_IDX = setdiff(1:4,[TOPPOINT_IDX LEFTPOINT_IDX RIGHTPOINT_IDX]);
            redPOINTS(2,:) = br(TOPPOINT_IDX,:);
            redPOINTS(4,:) = br(LEFTPOINT_IDX,:);
            redPOINTS(5,:) = br(MIDDLEPOINT_IDX,:);
            redPOINTS(6,:) = br(RIGHTPOINT_IDX,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            skeletonFrame = skeleton;
            skeletonFrame((redPOINTS(5,2) - 50):(redPOINTS(5,2) + 50),(redPOINTS(5,1) - 50):(redPOINTS(5,1) + 50)) = 0;
            skeletonFrame_lessTOP = skeletonFrame;
            skeletonFrame_lessTOP((redPOINTS(2,2) - 50):(redPOINTS(2,2) + 50),(redPOINTS(2,1) - 50):(redPOINTS(2,1) + 50)) = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % trace edge from left to top
            sk = [];
            [sk(:,2),sk(:,1)] = find(skeletonFrame);
            sourcePoint = snapTo(sk,redPOINTS(4,:));
            targetPoint = snapTo(sk,redPOINTS(2,:));
            ADJ = Radjacency(sk',2);
            [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
            path = sk(pathIDX,:);
            K = cwtK_filter(path,{7});
            K.K(1:50) = 0;
            K.K(end-50:end) = 0;
            [~,ptIDX] = max(abs(K.K));
            redPOINTS(1,:) = path(ptIDX,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % trace edge from top to right
            sk = [];
            [sk(:,2),sk(:,1)] = find(skeletonFrame);
            sourcePoint = snapTo(sk,redPOINTS(2,:));
            targetPoint = snapTo(sk,redPOINTS(6,:));
            ADJ = Radjacency(sk',2);
            [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
            path = sk(pathIDX,:);
            K = cwtK_filter(path,{7});
            K.K(1:50) = 0;
            K.K(end-50:end) = 0;
            [~,ptIDX] = min(K.K);
            redPOINTS(3,:) = path(ptIDX,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % trace edge from left to right
            sk = [];
            [sk(:,2),sk(:,1)] = find(skeletonFrame_lessTOP);
            sourcePoint = snapTo(sk,redPOINTS(4,:));
            targetPoint = snapTo(sk,redPOINTS(6,:));
            ADJ = Radjacency(sk',2);
            [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
            path = sk(pathIDX,:);
            K = cwtK_filter(path,{7});
            K.K(1:50) = 0;
            K.K(end-50:end) = 0;
            [~,ptIDX] = min(K.K);
            f1 = imdilate(K.K,ones(100,1)) == K.K;
            f2 = K.K > .02;
            fidx = find(f1.*f2);
            redPOINTS(7:8,:) = path(fidx,:);
            redPOINTS(7:8,:) = redPOINTS([8 7],:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure;
            imshow(qr,[]);
            hold on
            %plot(br(:,1),br(:,2),'bo');
            plot(redPOINTS(:,1),redPOINTS(:,2),'g*');
            %{
            CL = {'r*','m*','k*','c*'};
            for b = 1:size(boxCorners,3)
                for p = 1:size(boxCorners,1)
                    plot(boxCorners(p,1,b),boxCorners(p,2,b),CL{p})
                end
            end
            %}
            hold off;
            drawnow
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}


            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastRecord = size(QRTable,1) + 1;
            QRTable{lastRecord,'FileName'} = {qrFileName};
            QRTable{lastRecord,'cropBox'} = {{box}};
            for pt = 1:size(redPOINTS,1)
                QRTable{lastRecord,['red_point_' num2str(pt)]} = {redPOINTS(pt,:)};
            end
            %}

            imwrite(qr,qrFileName);
            csvwrite(strrep(qrFileName,'.tif','_box.csv'),box);
            csvwrite(strrep(qrFileName,'.tif','_redCorner.csv'),redPOINTS);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            deltaTM = toc;
            fprintf(['Ending QR data gather for:' num2str(e) ':' num2str(numel(FileList)) ':' num2str(mean(deltaTM)) '\n']);
            fprintf(['~remaining:' num2str(mean(deltaTM)*(numel(FileList) - e)/60/60/12) '...\n']);
        end
        
    catch ME
        getReport(ME)
    end

    
end
%% collect data for the red clicks
FilePath = '/mnt/tetra/nate/projectData/maizeSeedling/QR_mixed/';
FileList = {};
FileExt = {'tif'};
FileList = gdig(FilePath,FileList,FileExt,1);
FileList = FileList(randperm(numel(FileList)));
WZ = [21 21];
YS = [];
XS = [];
FR = 4;

for e = 1:500
    try
        tmpI = imread(FileList{e});
        tmpI = imresize(tmpI,1/FR);


        clickedPoints = readtext(strrep(FileList{e},'.tif','_redCorner.csv'));
        clickedPoints = cell2mat(clickedPoints);



        X = im2colF(double(tmpI),[WZ 3],[1 1 1]);
        Y = [];

        tmpX = reshape(X,[WZ 3 size(X,2)]);




        parfor pt = 1:8
            tmp = zeros(size(tmpI,1),size(tmpI,2));
            tmp(round(clickedPoints(pt,2)/FR),round(clickedPoints(pt,1)/FR)) = 1;
            tmp = imdilate(tmp,strel('disk',2,0));
            tmpY = im2colF(tmp,[WZ],[1 1]);
            Y(pt,:) = pt*tmpY((end-1)/2,:);
            %out = flattenMaskOverlay((tmpI+.1)/1.5,logical(tmp));
            %imshow(out,[]);
            %drawnow
            %pt
        end
        Y = sum(Y,1);

        fidx1 = find(Y~=0);
        fidx0 = find(Y==0);


        fidx0 = fidx0(randi(numel(fidx0),1,50));
        fidxT = [fidx0';fidx1'];

        YS = [YS Y(fidxT)];
        XS = [XS X(:,fidxT)];

        e
    catch ME
        getReport(ME)
    end
end
%% train red corner networks
X = reshape(XS,[WZ 3 size(XS,2)]);
%{
for e = 1:size(X,4)
    for k = 1:size(X,3)
        X(:,:,k,e) = bindVec(X(:,:,k,e));
    end
    e
    size(X,4)
end
%}
%
layers = [imageInputLayer([size(X,1) size(X,2) size(X,3)],'Normalization','none');
          convolution2dLayer(7,9);
          reluLayer();
          maxPooling2dLayer(3,'Stride',2);
          convolution2dLayer(3,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(9);
          softmaxLayer();
          classificationLayer()];
      
      
   
for e = 1
    fidx1 = find(YS~=0);
    fidx0 = find(YS==0);
    
    holdOutNumber1 = round(.80*numel(fidx1));
    holdOutNumber0 = round(.80*numel(fidx0));
    
    fidx1 = fidx1(randperm(numel(fidx1)));
    fidx0 = fidx0(randperm(numel(fidx0)));
    
    TRAIN1 = fidx1(1:holdOutNumber1);
    TEST1 = fidx1(holdOutNumber1:end);
    TRAIN0 = fidx0(1:holdOutNumber0);
    TEST0 = fidx0(holdOutNumber0:end);
    
    trainX = X(:,:,:,[TRAIN1 TRAIN0]);
    trainY = YS([TRAIN1 TRAIN0]);
    
    testX = X(:,:,:,[TEST1 TEST0]);
    testY = YS([TEST1 TEST0]);
   

    options = trainingOptions('sgdm','minibatch',512,'ValidationPatience',30,'Shuffle','every-epoch','ValidationFrequency',1000,'ValidationData',{testX,categorical(double(testY'))},'Plots','training-progress','MaxEpochs',100,'InitialLearnRate',.0001,'ExecutionEnvironment','parallel');

    redCornerNet = trainNetwork(trainX,categorical(double(trainY')),layers,options);
end
%% gather for conetainers

close all
FilePath = '/mnt/tetra/nate/projectData/maizeSeedling/NEF_to_tiff/';
FileList = {};
FileExt = {'tif'};
FileList = gdig(FilePath,FileList,FileExt,1);
FileList = FileList(randperm(numel(FileList)));
%% workbench
%%
%CS = [];
wsz = 17;
myF = zeros(2*wsz+1);
myF(((end-1)/2+1):end,:) = 1;
myS = myF;
myF(:,((end-1)/2+1):end) = 1;
myG = ones(2*wsz+1);
myF = [myS myF myG];

myG = [myG myG myG];
myF = [myF;myG];
myG = [myG;myG];
%rightPoints = {};
%leftPoints = {};
for e = 1:3000
            I = imread(FileList{e});
            I(:,(end-30):end,:) = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% start QR data gather
            Lab = rgb2lab(I);
            RED = Lab(:,:,2) > 25.5;
            RED = imclose(RED,strel('square',51));
            RED = bwlarge(RED);
            R = regionprops(RED);
            box = R(1).BoundingBox;
            box(1:2) = box(1:2) - pad;
            box(3:4) = box(3:4) + 2*pad;
            frame = imcrop(RED,box);
            frame = imclose(frame,strel('square',11));
            holes = imfill(frame,'holes') - frame;
            smallHoles = holes - bwareaopen(holes,1000);
            frame = logical(frame + smallHoles);
            qr = imcrop(I,box);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            oI = I;
            I(1:(box(2)+box(4)),:,:) = [];
            I((end-100):end,:,:) = [];
            if exist('GMModel_CLR')
                SZ = size(I);
                tI = imfilter(I,fspecial('gaussian',[21 21],[5]),'replicate');
                tI = rgb2lab(tI);
                tI = reshape(tI,[size(tI,1)*size(tI,2) size(tI,3)]);
                kidx = GMModel_CLR.cluster(double(tI));
                kidx = reshape(kidx,SZ(1:2));
                conemask = kidx==1;
                conemask = bwareaopen(conemask,3000);
                conemask = bwlarge(conemask,3);
                conemask = imfill(conemask,'holes');
                %imshow(conemask,[]);
                
                
                J = imfilter(double(conemask),myF,'replicate');
                K = imfilter(double(conemask),myG,'replicate');
                Q = (abs(K-J) < 2) & conemask ~= 0 & J >= .98*sum(myF(:));
                RL = regionprops(logical(Q),'Centroid');
                leftPoints{e} = [];
                for k = 1:numel(RL)
                    leftPoints{e} = [leftPoints{e};[RL(k).Centroid + [0 box(2)+box(4)]]];
                end
                
                
                J = imfilter(double(conemask),flipdim(myF,2),'replicate');
                K = imfilter(double(conemask),flipdim(myG,2),'replicate');
                QF = (abs(K-J) < 2) & conemask ~= 0 & J >= .98*sum(myF(:));
                RF = regionprops(logical(QF),'Centroid');
                rightPoints{e} = [];
                for k = 1:numel(RF)
                    rightPoints{e} = [rightPoints{e};[RF(k).Centroid + [0 box(2)+box(4)]]];
                end
                
                [q1 q2] = find(Q);
                [qf1 qf2] = find(QF);
                
              
                
                imshow(oI,[]);
                hold on
                if ~isempty(leftPoints{e}) & ~isempty(rightPoints{e})
                    
                    
                    leftPoints{e}(:,2) = leftPoints{e}(:,2) - 17;
                    rightPoints{e}(:,2) = rightPoints{e}(:,2) - 17;
                    
                    
                    plot(leftPoints{e}(:,1),leftPoints{e}(:,2),'g*');
                    plot(rightPoints{e}(:,1),rightPoints{e}(:,2),'r*');
                    hold on
                    
                    
                    
                    
                end
                %plot(q2,q1-wsz,'g*');
                %plot(qf2,qf1-wsz,'r*');
                drawnow
                hold off
                
            end
            I = rgb2lab(I);
            I = imresize(I,.2);
            I = reshape(I,[size(I,1)*size(I,2) size(I,3)]);
            I = I(randperm(size(I,1)),:);
            I = I(1:1000,:);
            CS = [CS;I];
            e
end
%% fit mixture model 
GMModel_CLR = fitgmdist(double(CS),3);
%% multiressample


BOX = [800 200];
BOXDET = [30 30];
BOXDET_normal = [800 800];
resSpace = linspace(.5,2,10);
resSpaceDET = 5;
%resSpaceDET = [1 2.5 5];
resSpace_normal = [.1 .5 1];
SZ = [100 400];
SZDET = [30 30];
SZDET_normal = [105 105];
X1_res = [];
X2_res = [];
Y_res1 = [];
Y_res2 = [];
OUTYR = [];
OUTYL = [];

X1_res_YES = [];
X1_res_NO = [];
X2_res_YES = [];
X2_res_NO = [];


X1_res_YES_normal = [];
X1_res_NO_normal = [];
X2_res_YES_normal = [];
X2_res_NO_normal = [];
applyNet = false;
gatherNet = true;


                
                REGION_LEFT_TRUE = [];
                REGION_LEFT_FALSE = [];
                REGION_RIGHT_TRUE = [];
                REGION_RIGHT_FALSE = [];
                
                 
                FINE_LEFT_TRUE = [];
                FINE_LEFT_FALSE = [];
                FINE_RIGHT_TRUE = [];
                FINE_RIGHT_FALSE = [];
               




[D] = generateMultiResSample(BOXDET_normal,resSpace_normal,SZDET_normal);
[regionDetectionDomain] = generateMultiResSampleDomain([800 800],[1],[20 20]);
[fineDetectionDomain] = generateMultiResSampleDomain([800 800],[.1 .5 1],[100 100]);
for e = 1:numel(rightPoints)
    
    if applyNet
        tic
        I = double(imread(FileList{e}));
        CROP_PAD = 200;
        G = rgb2gray(I/255);
        toView = G;
        for r = 1:4
            toView(1:CROP_PAD,:) = [];
            toView = imrotate(toView,90);
        end
        [s2 s1] = ndgrid(linspace(CROP_PAD+1,size(G,1)-CROP_PAD,NP(1)),linspace(CROP_PAD+1,size(G,2)-CROP_PAD,NP(2)));
        S = [s1(:) s2(:)];
        
        [toR] = sampleApply(G,S,regionDetectionDomain,[]);
        toL = coneDetNetL_region.predict(toR);     
        toL = reshape(toL(:,2),NP);
        toR = coneDetNetR_region.predict(toR);     
        toR = reshape(toR(:,2),NP);
        toR = imresize(toR,size(toView));
        toL = imresize(toL,size(toView));
        RRegion = toR > .98;
        LRegion = toL > .98;
        RRR = regionprops(RRegion,'PixelIdxList','Area');
        LLL = regionprops(LRegion,'PixelIdxList','Area');
        rightPT = [];
        leftPT = [];
        [rightPT(:,2),rightPT(:,1)] = find(RRegion);
        ridx = find(RRegion);
        lidx = find(LRegion);
        [leftPT(:,2),leftPT(:,1)] = find(LRegion);
        rightPT = bsxfun(@plus,rightPT,[400,400]);
        leftPT = bsxfun(@plus,leftPT,[400,400]);
        
        [toR_FINE] = sampleApply(G,rightPT,fineDetectionDomain,[]);
        [toL_FINE] = sampleApply(G,leftPT,fineDetectionDomain,[]);
        toRP = coneDetNetR_fine.predict(toR_FINE,'MiniBatchSize',2048);
        toLP = coneDetNetL_fine.predict(toL_FINE,'MiniBatchSize',2048);
        toR_FINE = zeros(size(toR));
        toR_FINE(ridx) = toRP(:,2);
        toL_FINE = zeros(size(toL));
        toL_FINE(lidx) = toLP(:,2);
        out = flattenMaskOverlay(toView,toR_FINE>.9,.7,'r');
        out = flattenMaskOverlay(out,toL_FINE>.9,.7,'b');
        
        [~,sidx1] = sort([LLL.Area],'descend');
        [~,sidx2] = sort([RRR.Area],'descend');
        LLL = LLL(sidx1(1:3));
        RRR = RRR(sidx2(1:3));
        for p = 1:3
            MSK = zeros(size(toView));
            MSK(LLL(p).PixelIdxList) = 1;
            MSK = MSK.*toL_FINE;
            MSK = logical(MSK > 0);
            tmp = regionprops(MSK,toL_FINE,'MeanIntensity','WeightedCentroid');
            [~,midx] = max(tmp.MeanIntensity);
            LP(p,:) = tmp(midx).WeightedCentroid;
        end
        for p = 1:3
            MSK = zeros(size(toView));
            MSK(RRR(p).PixelIdxList) = 1;
            MSK = MSK.*toR_FINE;
            MSK = logical(MSK > 0);
            tmp = regionprops(MSK,toL_FINE,'MeanIntensity','WeightedCentroid');
            [~,midx] = max(tmp.MeanIntensity);
            RP(p,:) = tmp(midx).WeightedCentroid;
        end
        
        out = imresize(out,.25);
        oSTORE{e} = out;
        toc
        
        
        
        I = imresize(I,1.25);
        Lab = rgb2lab(I);
        F = imfilter(Lab(:,:,1),fspecial('gaussian',[101 1],[11]),'replicate');
        [d1,d2] = gradient(F);
        sum_then_grad = sum(Lab(:,:,1),2);
        sum_then_grad = imfilter(sum_then_grad,fspecial('gaussian',[101 1],[31]),'replicate');
        sum_then_grad = gradient(sum_then_grad);
        grad_then_sum = sum(d2,2);
        grad_then_sum = imfilter(grad_then_sum,fspecial('gaussian',[101 1],[31]),'replicate');
        dv = mean([sum_then_grad grad_then_sum],2);
        [~,sidx] = min(dv);
        GR = rgb2gray(I);
        PAD_CONE = 300;
        TOPS = sidx-PAD_CONE;
        BOTS = min(sidx+PAD_CONE,size(I,1));
        SW = GR(TOPS:BOTS,:);
        
        
        GGR = GR;
        GGR(1:(3*end/4),:) = [];
        DS = 10;
        NP = [round(size(GGR,1)/DS) round(size(GGR,2)/DS)];
        Z = zeros(NP);
        J = GGR;
        for r = 1:4
            J(1:400,:) = [];
            J = imrotate(J,90);
        end
        domain = [];
        [s2 s1] = ndgrid(linspace(401,size(GGR,1)-400,NP(1)),linspace(401,size(GGR,2)-400,NP(2)));
        S = [s1(:) s2(:)];
        
        
        h1 = figure;
        h2 = figure;
        F = [];
        for s = 1:size(S,1)
            pt = S(s,:);
            
            %[f] = multiResSample(GGR,pt,resSpace_normal,BOXDET_normal,SZDET_normal);
            
            domain(:,1) = D(:,1) + pt(1);
            domain(:,2) = D(:,2) + pt(2);
            
            
            f = ba_interp2(double(GGR),domain(:,1),domain(:,2));
            f = reshape(f,[105 105 3]);
            F(:,:,:,s) = f;
            
            %tmp = coneDetNetR_normal.predict(f/255);
            %Z(s) = tmp(2);
            %{
            jZ = imresize(Z,size(J));
            imshow([(jZ) double(J)/255],[]);
            %}
            %drawnow
            
            %{
            figure(h1);
            imshow(Z,[]);
            figure(h2)
            imshow(GGR,[]);
            hold on
            %plot(domain(:,1),domain(:,2),'k.')
            plot(pt(1),pt(2),'r*')
            hold off
            drawnow
            %}
            s
            size(S,1)
        end
        SZF = size(F);
        rF = reshape(F,[prod(SZF(1:3)) SZF(4)]);
        cF = PCA_REPROJ_T(rF,qE,qU);
        cERR = sum((rF - PCA_BKPROJ_T(cF,qE,qU)).^2,1).^.5;
        ZZ = mvnpdf(cF',[0 0 0],qLAM);
        
        
        
        Z = coneDetNetR_normal.predict(F/255,'MiniBatchSize',2048);
        Z = reshape(Z(:,2),NP);
        jZ = imresize(Z,size(J));
        
        [j2,j1] = find(jZ > .5);
        [sidx] = find(jZ > .5);
        j1 = j1 + 400;
        j2 = j2 + 400;
        s1 = j1(mod(j1,5)==0 & mod(j2,5)==0);
        s2 = j2(mod(j1,5)==0 & mod(j2,5)==0);
        s1 = j1;
        s2 = j2;
        S = [s1(:) s2(:)];
        
        
        F = [];
        for s = 1:size(S,1)
            pt = S(s,:);
            
            %[f] = multiResSample(GGR,pt,resSpace_normal,BOXDET_normal,SZDET_normal);
            
            domain(:,1) = D(:,1) + pt(1);
            domain(:,2) = D(:,2) + pt(2);
            
            
            f = ba_interp2(double(GGR),domain(:,1),domain(:,2));
            f = reshape(f,[105 105 3]);
            F(:,:,:,s) = f;
            
            %tmp = coneDetNetR_normal.predict(f/255);
            %Z(s) = tmp(2);
            %{
            jZ = imresize(Z,size(J));
            imshow([(jZ) double(J)/255],[]);
            %}
            %drawnow
            
            %{
            figure(h1);
            imshow(Z,[]);
            figure(h2)
            imshow(GGR,[]);
            hold on
            %plot(domain(:,1),domain(:,2),'k.')
            plot(pt(1),pt(2),'r*')
            hold off
            drawnow
            %}
            s
            size(S,1)
        end
        
        Z = coneDetNetR_normal.predict(F/255,'MiniBatchSize',2048);
        jZ2 = zeros(size(jZ));
        jZ2(sidx) = Z(:,2);
        
        SW = imresize(SW,resSpaceDET^-1);
        SW_SL = im2colF(double(SW),BOXDET,[1 1]);
        SW_SL = reshape(SW_SL,[BOXDET 1 size(SW_SL,2)]);
        QRDL = coneDetNetL.predict(SW_SL/255,'MiniBatchSize',2048);
        QRDR = coneDetNetR.predict(SW_SL/255,'MiniBatchSize',2048);
        HOPEL = col2im(QRDL(:,2),BOXDET,size(SW),'sliding'); 
        HOPER = col2im(QRDR(:,2),BOXDET,size(SW),'sliding');
        HOPE = cat(3,HOPEL,zeros(size(HOPER)),HOPER);
        
        
        
        
        RR = regionprops(HOPER>.95,'Centroid');
        toSample = rgb2gray(I);
        BUF = 100;
        LOOP = 5;
        PT = zeros(LOOP,2,numel(RR));
        isV = zeros(LOOP,numel(RR));
        for re = 1:numel(RR)
           
            PT(1,:,re) = resSpaceDET*(RR(re).Centroid + BOXDET/2) +[0 TOPS];
            resPredict(1,re) = 1;
            for loop = 1:LOOP
                [dd2,dd1] = ndgrid((PT(loop,2,re)-BUF):4:(PT(loop,2,re)+BUF),(PT(loop,1,re)-BUF):4:(PT(loop,1,re)+BUF));
                dd = [dd1(:),dd2(:)];
                tmpO = [];
                for pt = 1:numel(dd1)
                    tmpBOX = BOXDET_normal;%*resPredict(loop,re);
                    [tmpy] = multiResSample(toSample,dd(pt,:),resSpace_normal,tmpBOX,SZDET_normal);
                    tmpO(:,:,:,pt) = squeeze(tmpy);
                    %{
                    Q_normal_test = coneDetNetR_normal.predict(tmpO(:,:,:,pt)/255);
                    imshow(tmpO(:,:,:,pt)/255,[]);
                    hold on
                    plot(250/2,250/2,'r*');
                    title(num2str(Q_normal_test))
                    hold off
                    drawnow
                    %}
                end
                Q_normal = coneDetNetR_normal.predict(tmpO/255,'MiniBatchSize',2048);
                Q_normal = reshape(Q_normal(:,2),size(dd1));
                [isV(loop+1),pidx] = max(Q_normal(:));
                pidx = find(Q_normal(:) == isV(loop+1));
                PT(loop+1,:,re) = mean(dd(pidx,:),1);
                
                
                
                
                RES_SAM = multiResSample(toSample,PT(loop+1,:,re),1,BOX,SZ);
                resPredict(loop+1,re) = coneResNet.predict(RES_SAM/255) + Y_res1_U;
            end
            
            
            
            imshow(toSample,[]);
            hold on
            plot(PT(:,1,re),PT(:,2,re),'b')
            plot(PT(1,1,re),PT(1,2,re),'g*')
            plot(PT(end,1,re),PT(end,2,re),'r*')
            waitforbuttonpress
           
        end
        
        
        
        
        
    end
    
    
    if gatherNet
        if true%contains(FileList{e},'MN')
            if size(rightPoints{e},1) == 3 && size(leftPoints{e},1)  == 3
                OUTR = [];
                OUTL = [];
                OUTYR = [];
                OUTYL = [];
                OUTR_YES = [];
                OUTR_NO = [];
                OUTL_YES = [];
                OUTL_NO = [];
                OUTR_YES_normal = [];
                OUTR_NO_normal = [];
                OUTL_YES_normal = [];
                OUTL_NO_normal = [];
                
                I = imread(FileList{e});
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample for classifiction - low res
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                trueMask = zeros(size(I,1),size(I,2));
                for pt = 1:size(rightPoints{e},1)
                    trueMask(round(rightPoints{e}(pt,2)),round(rightPoints{e}(pt,1))) = 1;
                end
                [TrightRegionPoints,FrightRegionPoints] = generateRegionPoints(trueMask,20,.25,10,-1,500);
                
                trueMask = zeros(size(I,1),size(I,2));
                for pt = 1:size(leftPoints{e},1)
                    trueMask(round(leftPoints{e}(pt,2)),round(leftPoints{e}(pt,1))) = 1;
                end
                
                [TleftRegionPoints,FleftRegionPoints] = generateRegionPoints(trueMask,20,.25,10,-1,500);
                
                [TsampleLeftRegions] = sampleApply(double(rgb2gray(I)),TleftRegionPoints,regionDetectionDomain,[]);
                [TsampleRightRegions] = sampleApply(double(rgb2gray(I)),TrightRegionPoints,regionDetectionDomain,[]);
                [FsampleLeftRegions] = sampleApply(double(rgb2gray(I)),FleftRegionPoints,regionDetectionDomain,[]);
                [FsampleRightRegions] = sampleApply(double(rgb2gray(I)),FrightRegionPoints,regionDetectionDomain,[]);
                
                REGION_LEFT_FALSE = cat(4,REGION_LEFT_FALSE,FsampleLeftRegions);
                REGION_LEFT_TRUE = cat(4,REGION_LEFT_TRUE,TsampleLeftRegions);
                REGION_RIGHT_FALSE = cat(4,REGION_RIGHT_FALSE,FsampleRightRegions);
                REGION_RIGHT_TRUE = cat(4,REGION_RIGHT_TRUE,TsampleRightRegions);
                
                
                
                
                
                
                trueMask = zeros(size(I,1),size(I,2));
                for pt = 1:size(rightPoints{e},1)
                    trueMask(round(rightPoints{e}(pt,2)),round(rightPoints{e}(pt,1))) = 1;
                end
                [TrightFinePoints,FrightFinePoints] = generateRegionPoints(trueMask,5,1,1,-1,100);
                
                trueMask = zeros(size(I,1),size(I,2));
                for pt = 1:size(leftPoints{e},1)
                    trueMask(round(leftPoints{e}(pt,2)),round(leftPoints{e}(pt,1))) = 1;
                end
                
                [TleftFinePoints,FleftFinePoints] = generateRegionPoints(trueMask,5,1,1,-1,100);
                
                [TsampleLeftFine] = sampleApply(double(rgb2gray(I)),TleftFinePoints,fineDetectionDomain,[]);
                [TsampleRightFine] = sampleApply(double(rgb2gray(I)),TrightFinePoints,fineDetectionDomain,[]);
                [FsampleLeftFine] = sampleApply(double(rgb2gray(I)),FleftFinePoints,fineDetectionDomain,[]);
                [FsampleRightFine] = sampleApply(double(rgb2gray(I)),FrightFinePoints,fineDetectionDomain,[]);
                
                FINE_LEFT_FALSE = cat(4,FINE_LEFT_FALSE,FsampleLeftFine);
                FINE_LEFT_TRUE = cat(4,FINE_LEFT_TRUE,TsampleLeftFine);
                FINE_RIGHT_FALSE = cat(4,FINE_RIGHT_FALSE,FsampleRightFine);
                FINE_RIGHT_TRUE = cat(4,FINE_RIGHT_TRUE,TsampleRightFine);
                
                
                
                
                
                %{
                
                for pt = 1:size(rightPoints{e},1)
                    PT = (rightPoints{e}(pt,:));
                    toSample = rgb2gray(I);
                    [tmpO_R_BINARY_YES] = multiResSample(toSample,PT,resSpaceDET,BOXDET,SZDET);
                    tmpO_R_BINARY_YES = squeeze(tmpO_R_BINARY_YES);
                    try
                        for sam = 1:10
                            NOPT = PT + 400*(rand(1,2)+.1);
                            [tmpO_R_BINARY_NO] = multiResSample(toSample,NOPT,resSpaceDET,BOXDET,SZDET);
                            tmpO_R_BINARY_NO = squeeze(tmpO_R_BINARY_NO);
                            OUTR_NO = cat(4,OUTR_NO,tmpO_R_BINARY_NO);
                        end
                    catch
                    end
                    OUTR_YES = cat(4,OUTR_YES,tmpO_R_BINARY_YES);
                    
                end
                for pt = 1:size(leftPoints{e},1)
                    PT = (leftPoints{e}(pt,:));
                    toSample = rgb2gray(I);
                    [tmpO_L_BINARY_YES] = multiResSample(toSample,PT,resSpaceDET,BOXDET,SZDET);
                    tmpO_L_BINARY_YES = squeeze(tmpO_L_BINARY_YES);
                    try
                        for sam = 1:10
                            NOPT = PT + 400*(rand(1,2)+.1);
                            [tmpO_L_BINARY_NO] = multiResSample(toSample,NOPT,resSpaceDET,BOXDET,SZDET);
                            tmpO_L_BINARY_NO = squeeze(tmpO_L_BINARY_NO);
                            OUTL_NO = cat(4,OUTL_NO,tmpO_L_BINARY_NO);
                        end
                    catch
                    end

                    OUTL_YES = cat(4,OUTL_YES,tmpO_L_BINARY_YES);
                   
                end
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample for classifiction
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %{
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample for classifiction - normal res
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for pt = 1:size(rightPoints{e},1)
                    PT = (rightPoints{e}(pt,:));
                    toSample = rgb2gray(I);
                    [tmpO_R_BINARY_YES_normal] = multiResSample(toSample,PT,resSpace_normal,BOXDET_normal,SZDET_normal);
                    tmpO_R_BINARY_YES_normal = squeeze(tmpO_R_BINARY_YES_normal);
                    try
                        for sam = 1:20
                            NOPT = PT + 400*(rand(1,2)+.1);
                            [tmpO_R_BINARY_NO_normal] = multiResSample(toSample,NOPT,resSpace_normal,BOXDET_normal,SZDET_normal);
                            tmpO_R_BINARY_NO_normal = squeeze(tmpO_R_BINARY_NO_normal);
                            OUTR_NO_normal = cat(4,OUTR_NO_normal,tmpO_R_BINARY_NO_normal);
                        end
                    catch
                    end
                    
                    try
                        for sam = 1:20
                            NOPT = PT + 10*(rand(1,2)+.1);
                            [tmpO_R_BINARY_NO_normal] = multiResSample(toSample,NOPT,resSpace_normal,BOXDET_normal,SZDET_normal);
                            tmpO_R_BINARY_NO_normal = squeeze(tmpO_R_BINARY_NO_normal);
                            OUTR_NO_normal = cat(4,OUTR_NO_normal,tmpO_R_BINARY_NO_normal);
                        end
                    catch
                    end
                    OUTR_YES_normal = cat(4,OUTR_YES_normal,tmpO_R_BINARY_YES_normal);
                    
                end
                for pt = 1:size(leftPoints{e},1)
                    PT = (leftPoints{e}(pt,:));
                    toSample = rgb2gray(I);
                    [tmpO_L_BINARY_YES_normal] = multiResSample(toSample,PT,resSpace_normal,BOXDET_normal,SZDET_normal);
                    tmpO_L_BINARY_YES_normal = squeeze(tmpO_L_BINARY_YES_normal);
                    try
                        for sam = 1:20
                            NOPT = PT + 400*(rand(1,2)+.1);
                            [tmpO_L_BINARY_NO_normal] = multiResSample(toSample,NOPT,resSpace_normal,BOXDET_normal,SZDET_normal);
                            tmpO_L_BINARY_NO_normal = squeeze(tmpO_L_BINARY_NO_normal);
                            OUTL_NO_normal = cat(4,OUTL_NO_normal,tmpO_L_BINARY_NO_normal);
                        end
                    catch
                    end
                    
                    try
                        for sam = 1:20
                            NOPT = PT + 10*(rand(1,2)+.1);
                            [tmpO_L_BINARY_NO_normal] = multiResSample(toSample,NOPT,resSpace_normal,BOXDET_normal,SZDET_normal);
                            tmpO_L_BINARY_NO_normal = squeeze(tmpO_L_BINARY_NO_normal);
                            OUTL_NO_normal = cat(4,OUTL_NO_normal,tmpO_L_BINARY_NO_normal);
                        end
                    catch
                    end
                    

                    OUTL_YES_normal = cat(4,OUTL_YES_normal,tmpO_L_BINARY_YES_normal);
                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample for classifiction
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}

                
                %{
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample for resolution
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for pt = 1:size(rightPoints{e},1)
                    PT = (rightPoints{e}(pt,:));
                    toSample = rgb2gray(I);
                    [tmpO] = multiResSample(toSample,PT,resSpace,BOX,SZ);
                    %{
                    if exist('coneResNet')
                        resPreR = coneResNet.predict(tmpO/255) + Y_res1_U;
                    end
                    %}
                    OUTR = cat(4,OUTR,tmpO);
                    OUTYR = [OUTYR;resSpace'];
                end
                for pt = 1:size(leftPoints{e},1)
                    PT = (leftPoints{e}(pt,:));
                    toSample = rgb2gray(I);
                    [tmpO] = multiResSample(toSample,PT,resSpace,BOX,SZ);
                    OUTL = cat(4,OUTL,tmpO);
                    OUTYL = [OUTYL;resSpace'];
                end
                X1_res = cat(4,X1_res,OUTR);
                X2_res = cat(4,X2_res,OUTL);


                X1_res_YES = cat(4,X1_res_YES,OUTR_YES);
                X1_res_NO = cat(4,X1_res_NO,OUTR_NO);

                X2_res_YES = cat(4,X2_res_YES,OUTL_YES);
                X2_res_NO = cat(4,X2_res_NO,OUTL_NO);
                
                
                
                X1_res_YES_normal = cat(4,X1_res_YES_normal,OUTR_YES_normal);
                X1_res_NO_normal = cat(4,X1_res_NO_normal,OUTR_NO_normal);

                X2_res_YES_normal = cat(4,X2_res_YES_normal,OUTL_YES_normal);
                X2_res_NO_normal = cat(4,X2_res_NO_normal,OUTL_NO_normal);

                Y_res1 = [Y_res1;OUTYR];
                Y_res2 = [Y_res2;OUTYL];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample for resolution
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %}
                
                
                
            end
        end
         e
    end
end
%%
cropPackage = {};
for e = 1:500
    tic
    cropPackage{e} = getSeedlingCropBoxes(FileList{9},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine);
    toc
end

%% train region
SZ1 = size(REGION_LEFT_TRUE);
%SZ1 = size(X_CL1_normal_n);
layers = [ ...
    imageInputLayer([SZ1(1:3)],'Normalization','none')
    convolution2dLayer(5,5)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
options = trainingOptions('sgdm','Shuffle','every-epoch','Plots','training-progress','MaxEpochs',10,'InitialLearnRate',0.001,'ExecutionEnvironment','parallel');
X_CL1_normal = cat(4,REGION_LEFT_TRUE,REGION_LEFT_FALSE);
Y_CL1_normal = [ones(size(REGION_LEFT_TRUE,4),1);0*ones(size(REGION_LEFT_FALSE,4),1)];
coneDetNetL_region = trainNetwork(X_CL1_normal/255,categorical(Y_CL1_normal),layers,options);
X_CL1_normal = cat(4,REGION_RIGHT_TRUE,REGION_RIGHT_FALSE);
Y_CL1_normal = [ones(size(REGION_RIGHT_TRUE,4),1);0*ones(size(REGION_RIGHT_FALSE,4),1)];
coneDetNetR_region = trainNetwork(X_CL1_normal/255,categorical(Y_CL1_normal),layers,options);
%% train fine
SZ1 = size(FINE_LEFT_TRUE);
%SZ1 = size(X_CL1_normal_n);
layers = [ ...
    imageInputLayer([SZ1(1:3)],'Normalization','none')
    convolution2dLayer(13,3)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
options = trainingOptions('sgdm','Shuffle','every-epoch','Plots','training-progress','MaxEpochs',10,'InitialLearnRate',0.001,'ExecutionEnvironment','parallel');
X_CL1_normal = cat(4,FINE_LEFT_TRUE,FINE_LEFT_FALSE);
Y_CL1_normal = [ones(size(FINE_LEFT_TRUE,4),1);0*ones(size(FINE_LEFT_FALSE,4),1)];
X_CL1_normal = X_CL1_normal/255;
coneDetNetL_fine = trainNetwork(X_CL1_normal,categorical(Y_CL1_normal),layers,options);
X_CL1_normal = cat(4,FINE_RIGHT_TRUE,FINE_RIGHT_FALSE);
Y_CL1_normal = [ones(size(FINE_RIGHT_TRUE,4),1);0*ones(size(FINE_RIGHT_FALSE,4),1)];
coneDetNetR_fine = trainNetwork(X_CL1_normal/255,categorical(Y_CL1_normal),layers,options);

%% train classification layer
SZ1 = size(X1_res_YES_normal);
%SZ1 = size(X_CL1_normal_n);
layers = [ ...
    imageInputLayer([SZ1(1:3)],'Normalization','none')
    convolution2dLayer(13,5)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
options = trainingOptions('sgdm','Shuffle','every-epoch','Plots','training-progress','MaxEpochs',10,'InitialLearnRate',0.001,'ExecutionEnvironment','parallel');
X_CL1_normal = cat(4,X1_res_YES_normal,X1_res_NO_normal,X2_res_YES_normal);
Y_CL1_normal = [ones(size(X1_res_YES_normal,4),1);0*ones(size(X1_res_NO_normal,4),1);0*ones(size(X2_res_YES_normal,4),1)];

%X_CL1_normal = cat(4,X1_res_YES_normal,X1_res_NO_normal);
%Y_CL1_normal = [ones(size(X1_res_YES_normal,4),1);0*ones(size(X1_res_NO_normal,4),1)];


%{
    SZQ = size(X1_res_YES_normal);
    R = reshape(X1_res_YES_normal,[prod(SZQ(1:3)) SZQ(4)]);
    [qS qC qU qE qL qERR qLAM] = PCA_FIT_FULL_T(R,3);

%}

coneDetNetR_normal = trainNetwork(X_CL1_normal/255,categorical(Y_CL1_normal),layers,options);

%{

    for e = 1:size(X_CL1_normal,4)
        X_CL1_normal_n(:,:,:,:,e) = imresize(X_CL1_normal(:,:,:,e),[55 55]);
    end
%}



X_CL2_normal = cat(4,X2_res_YES_normal,X2_res_NO_normal,X1_res_YES_normal);
Y_CL2_normal = [ones(size(X2_res_YES_normal,4),1);0*ones(size(X2_res_NO_normal,4),1);0*ones(size(X1_res_YES_normal,4),1)];

%X_CL2_normal = cat(4,X2_res_YES_normal,X2_res_NO_normal);
%Y_CL2_normal = [ones(size(X2_res_YES_normal,4),1);0*ones(size(X2_res_NO_normal,4),1)];

coneDetNetL_normal = trainNetwork(X_CL2_normal/255,categorical(Y_CL2_normal),layers,options);


SZ1 = size(X1_res_YES);
layers = [ ...
    imageInputLayer([SZ1(1:3)],'Normalization','none')
    convolution2dLayer(13,3)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
options = trainingOptions('sgdm','Shuffle','every-epoch','Plots','training-progress','MaxEpochs',500,'InitialLearnRate',0.001,'ExecutionEnvironment','parallel');
X_CL1 = cat(4,X1_res_YES,X1_res_NO,X2_res_YES);
Y_CL1 = [ones(size(X1_res_YES,4),1);0*ones(size(X1_res_NO,4),1);0*ones(size(X2_res_YES,4),1)];
coneDetNetR = trainNetwork(X_CL1/255,categorical(Y_CL1),layers,options);
X_CL2 = cat(4,X2_res_YES,X2_res_NO,X1_res_YES);
Y_CL2 = [ones(size(X2_res_YES,4),1);0*ones(size(X2_res_NO,4),1);0*ones(size(X1_res_YES,4),1)];
coneDetNetL = trainNetwork(X_CL2/255,categorical(Y_CL2),layers,options);
% train classification layer

%% train res networks
SZ1 = size(X1_res);
layers = [ ...
    imageInputLayer([SZ1(1:3)],'Normalization','none')
    convolution2dLayer(12,10)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','Shuffle','every-epoch','Plots','training-progress','MaxEpochs',10,'InitialLearnRate',0.001,'ExecutionEnvironment','parallel');
Y_res1_U = mean(Y_res1);
coneResNet = trainNetwork(X1_res/255,Y_res1 - Y_res1_U,layers,options);
%% try resolution networks
close all
y_pre1 = coneResNet.predict(X1_res/255);
plot(y_pre1+Y_res1_U,Y_res1,'.')
%%
[conePoint coneUse] = zoomGather_nonDB(FileList(1:50),[200 200],true);
%%
close all
points = {};
for e = 1:3
    [~,BOX{e}] = imcrop(I);
end
%%
for e = 1:3
    points{e} = detectMinEigenFeatures(rgb2gray(I),'ROI',BOX{e});
    points{e} = selectStrongest(points{e},29);
    [features{e},validPoints{e}] = extractFeatures(rgb2gray(I),points{e},'Method','Block','BlockSize',17);
    
    %[features{e},validPoints{e}] = extractHOGFeatures(I,points{e});
    
    
    
    v = nchoosek(1:size(features{e},1),2);
    delta = pdist(validPoints{e}.Location);
    Sfeatures{e} = [[features{e}(v(:,1),:) features{e}(v(:,2),:) delta'];[features{e}(v(:,2),:) features{e}(v(:,1),:) delta']];
    PP{e} = [[validPoints{e}.Location(v(:,1),:) validPoints{e}.Location(v(:,2),:)];[validPoints{e}.Location(v(:,2),:) validPoints{e}.Location(v(:,1),:)]];






end

imshow(I,[]);
hold all
for e = 1:numel(points)
  plot(validPoints{e}.selectStrongest(50));
  hold all
end
waitforbuttonpress
v = nchoosek(1:3,2);
for e = 1:size(v,1)
    [indexPairs{e},matchmetric{e}] = matchFeatures(features{v(e,1)},features{v(e,2)},'Unique',true,'MatchThreshold',10);
    for p = 1:size(indexPairs{e},1)
        seg = [validPoints{v(e,1)}.Location(indexPairs{e}(p,1),:);validPoints{v(e,2)}.Location(indexPairs{e}(p,2),:)];
        imshow(I,[]);
        hold on
        plot(seg(:,1),seg(:,2),'c');
        hold off
        waitforbuttonpress
    end
end
%%

v = nchoosek(1:3,2);
for e = 1:size(v,1)
    
    DIS = pdist2(Sfeatures{v(e,1)},Sfeatures{v(e,2)});
    [~,IDX2] = min(DIS,[],2);
    [~,IDX1] = min(DIS,[],1);
    IDX = find(IDX1'==IDX2);
    indexPairs{e} = IDX';
    
    [indexPairs{e},matchmetric{e}] = matchFeatures(Sfeatures{v(e,1)},Sfeatures{v(e,2)},'Unique',true,'MatchThreshold',1,'Metric','SAD','MaxRatio',.4);
    
    
    for p = 1:size(indexPairs{e},1)
        seg1 = reshape(PP{v(e,1)}(indexPairs{e}(p,1),:),[2 2])';
        seg2 = reshape(PP{v(e,2)}(indexPairs{e}(p,2),:),[2 2])';
        imshow(I,[]);
        hold on
        plot(seg1(:,1),seg1(:,2),'c');
        plot(seg2(:,1),seg2(:,2),'c');
        hold off
        waitforbuttonpress
    end
end
%%
bottomSTACK = [];
conetainerBUMP_NEG = 350;
conetainerBUMP = 55;
HMC = {};
for e = 1:3000


    tic
    I = double(imread(FileList{e}))/255;
    
    
    if e == 1
        WIDTH = size(I,2);
    end
    I(:,(end-100):end,:) = [];




    %%%% start QR data gather
    Lab = rgb2lab(I);
    RED = Lab(:,:,2) > 25.5;
    RED = imclose(RED,strel('square',51));
    RED = bwlarge(RED);
    R = regionprops(RED);
    box = R(1).BoundingBox;
    box(1:2) = box(1:2) - pad;
    box(3:4) = box(3:4) + 2*pad;
     %%%% start QR data gather



    %%%% gather for conetainers
    Lab = rgb2lab(I);
    F = imfilter(Lab(:,:,1),fspecial('gaussian',[101 1],[11]),'replicate');
    [d1,d2] = gradient(F);
    sum_then_grad = sum(Lab(:,:,1),2);
    sum_then_grad = imfilter(sum_then_grad,fspecial('gaussian',[101 1],[31]),'replicate');
    sum_then_grad = gradient(sum_then_grad);
    grad_then_sum = sum(d2,2);
    grad_then_sum = imfilter(grad_then_sum,fspecial('gaussian',[101 1],[31]),'replicate');
    dv = mean([sum_then_grad grad_then_sum],2);
    [~,sidx] = min(dv);
    sidx = sidx - conetainerBUMP;
    sidx2 = sidx + conetainerBUMP_NEG;
    sidx2 = min(sidx2,size(I,1));
    
    close all
    figure;
    imshow(I,[]);
    hold on
    plot(1:size(I,2),sidx*ones(1,size(I,2)),'r');
    plot(1:size(I,2),sidx2*ones(1,size(I,2)),'r');
    drawnow
    hold off
    
    
    bottom = I(sidx:sidx2,:,:);
    bottom = imresize(bottom,[100 WIDTH]);
    bottom = im2colF(rgb2gray(bottom),[size(bottom,1) 50],[1 1]);
    
    
    if exist('GMModel_CONE')
       
        %gbottom = rgb2gray(bottom);
        pp = PCA_REPROJ_T(bottom,bE,bU);
        pp = padarray(pp,[0 25],'both','replicate');
        ppDISPLAY = interp1(linspace(1,size(pp,2),size(pp,2))',pp',linspace(1,size(pp,2),size(I,2))','nearest');
        ppDISPLAY = ppDISPLAY';
        hold on
        CL = {'r' 'g' 'b'};
        for k = 1:size(ppDISPLAY,1)
            plot(1:size(I,2),ppDISPLAY(k,:)*10 + size(I,1)*.25,CL{k})
        end
        
        
        pp = GMModel_CONE.cluster(pp');
        pp = interp1(linspace(1,size(pp,1),size(pp,1))',pp,linspace(1,size(pp,1),size(I,2))','nearest');
        
        
        
        
       
        plot(1:size(I,2),pp*100 + size(I,1)*.75,'m')
        hold off
        drawnow
    end
    
    if exist('hmm')
        %gbottom = rgb2gray(bottom);
        pp = PCA_REPROJ_T(bottom,bE,bU);
        [states,prob] = hmm.Viterbi(pp,[1 1]);
        states = padarray(states,[0 25],'both','replicate');
      
        
        
        HMC{e} = [states;pp];
        
        
        states = interp1(linspace(1,size(states',1),size(states',1))',states',linspace(1,size(states',1),size(I,2))','nearest');
       
        hold on
        plot(1:size(I,2),states*100 + size(I,1)*.5,'g')
        hold off
        drawnow
    end
    
    
    
    
    bottomSTACK = cat(4,bottomSTACK,bottom);
    e
end
%% rgb 2 gray
%{
 bottomSTACK_gray = [];
for e = 1:size(bottomSTACK,4)
    bottomSTACK_gray(:,:,e) = rgb2gray(bottomSTACK(:,:,:,e));
    e
    size(bottomSTACK,4)
end
%}
bottomSTACK_gray = squeeze(bottomSTACK);
%%
SZ = size(bottomSTACK_gray);
bottomSTACK_gray_R = reshape(bottomSTACK_gray,[SZ(1) prod(SZ(2:3))]);
[~,bC,bU,bE] = PCA_FIT_FULL_T(bottomSTACK_gray_R,2);
options = statset('Display','iter');
GMModel_CONE = fitgmdist([bC'],3,'Options',options);
%%

%% finding conetainer - step 2 make hmm 
kidx = GMModel_CONE.cluster(bC');
backGroundIDX = find(kidx==1);
coneIDX = find(kidx==2);
smallRE = regionprops(kidx==2,'Area','PixelIdxList');
smallRE([smallRE.Area] > 250) = [];
for e = 1:numel(smallRE)
    kidx(smallRE(e).PixelIdxList) = 3;
end
coneIDX = find(kidx==2);
containerIDX = find(kidx==3);

kidx(coneIDX(bC(1,coneIDX) > -10)) = 3;



coneIDX = find(kidx==2);
containerIDX = find(kidx==3);

background_node = hmm_node('background');
cone_node = hmm_node('cone');
container_node_PRE = hmm_node('container_pre');
container_node_POST = hmm_node('container_post');

backgroundCov = cov(bC(:,backGroundIDX)');
backgroundMean = mean(bC(:,backGroundIDX)');

coneCov = cov(bC(:,coneIDX)');
coneMean = mean(bC(:,coneIDX)');

conetainerCov = cov(bC(:,containerIDX)');
conetainerMean = mean(bC(:,containerIDX)');


backgroundDis = myProb(backgroundMean,backgroundCov);
coneDis = myProb(coneMean,coneCov);
conetainerDis = myProb(conetainerMean,conetainerCov);
    
background_node.attachDistribution(backgroundDis,1);
cone_node.attachDistribution(coneDis,1);
container_node_PRE.attachDistribution(conetainerDis,1);
container_node_POST.attachDistribution(conetainerDis,1);

backGroundTobackGround = constantTransitionFunction(.5);
backGroundTobackGround = VarheavisideTransitionFunction(200,@(x,y)lt(x,y),.9,1);
backGroundTocontainer = constantTransitionFunction(.001);
backGroundTocontainer = VarheavisideTransitionFunction(200,@(x,y)ge(x,y),.9,1);

containerTocontainer = constantTransitionFunction(.8);
containerTocontainer = VarheavisideTransitionFunction(400,@(x,y)lt(x,y),.9,1);
%containerTocontainer = heavisideTransitionFunction(550,@(x,y)lt(x,y));
%containerTocontainer = heavisideTransitionFunction(300,@(x,y)lt(x,y));
%containerPRETocontainerPOST = constantTransitionFunction(.0001);
containerTobackground = constantTransitionFunction(.1);
%containerTobackground = VarheavisideTransitionFunction(400,@(x,y)lt(x,y),.9999,1);
%containerTobackground = heavisideTransitionFunction(550,@(x,y)ge(x,y));
%containerTobackground = heavisideTransitionFunction(300,@(x,y)ge(x,y));
%containerTocone = constantTransitionFunction(.00000001);
containerTocone = VarheavisideTransitionFunction(400,@(x,y)lt(x,y),.01,1);
%containerTocone = heavisideTransitionFunction(550,@(x,y)ge(x,y));

coneTocone = constantTransitionFunction(.9);
coneTocone = heavisideTransitionFunction(215,@(x,y)lt(x,y));
coneTocone = VarheavisideTransitionFunction(300,@(x,y)lt(x,y),1,1);
coneToconetainer = constantTransitionFunction(.001);
coneToconetainer = VarheavisideTransitionFunction(300,@(x,y)ge(x,y),1,1);
%coneToconetainer = heavisideTransitionFunction(150,@(x,y)ge(x,y));

background_node.attachNode(background_node,backGroundTobackGround);
background_node.attachNode(container_node_PRE,backGroundTocontainer);

container_node_PRE.attachNode(container_node_PRE,containerTocontainer);
container_node_PRE.attachNode(cone_node,containerTocone);
%container_node_PRE.attachNode(container_node_POST,containerPRETocontainerPOST);

container_node_POST.attachNode(container_node_POST,containerTocontainer);
container_node_POST.attachNode(background_node,containerTobackground);

cone_node.attachNode(cone_node,coneTocone);
cone_node.attachNode(container_node_POST,coneToconetainer);


hmm = my_hmm();
hmm.addNode(background_node);
hmm.addNode(container_node_PRE);
hmm.addNode(cone_node);
hmm.addNode(container_node_POST);

hmm.dn = ones(size(I,2),1);
%%
hmm.update(HMC,[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auto clicks for RED CORNER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
autoPath = '/mnt/tetra/nate/autoRCNN/';
close all
pad = 50;
QRSHEET = [];
QR_PIECES = [];

cnt = 1;
clickedFileList = {};
clickPoints = {};
preD = false;
bottomSTACK = [];
TOP_STOCK = [];
conetainerBUMP = 50;
conetainerBUMP_NEG = 400;
HMC = {};
MASTER_PLANT_STACK = [];
%FileList = FileList(randperm(numel(FileList)));
gatherQR = false;
CROP_BUMP = 150;
for e = 1:1:200%numel(FileList)
    try
        
       
        
        tic
        I = double(imread(FileList{e}))/255;
        I(:,(end-100):end,:) = [];
        close all

        
        
         
         %%%% start QR data gather
            Lab = rgb2lab(I);
            RED = Lab(:,:,2) > 25.5;
            RED = imclose(RED,strel('square',51));
            RED = bwlarge(RED);
            R = regionprops(RED);
            box = R(1).BoundingBox;
            box(1:2) = box(1:2) - pad;
            box(3:4) = box(3:4) + 2*pad;
         %%%% start QR data gather
        
        
        
        %%%% gather for conetainers
        Lab = rgb2lab(I);
        vec2 = sum(Lab(:,:,1),2);
        svec2 = imfilter(vec2,fspecial('gaussian',[101 1],[21]),'replicate');
        dv = gradient(svec2);
        [~,sidx ] = min(dv);
        sidx = sidx - conetainerBUMP;
        sidx2 = sidx + conetainerBUMP_NEG;
        figure;
        imshow(I,[]);
        hold on
        plot(1:size(I,2),sidx*ones(1,size(I,2)),'r');
        plot(1:size(I,2),sidx2*ones(1,size(I,2)),'r');
        
        bottom = I(sidx:sidx2,:,:);
        bottom = imresize(bottom,[100 size(I,2)]);
        bottom = bsxfun(@minus,rgb2gray(bottom),mean(rgb2gray(bottom),2));
        bottomSTACK = [bottomSTACK bottom];
        
        %bottom = padarray(bottom,[0 300],'both','replicate');
        
        
        vecB = PCA_REPROJ_T((bottom),bE,bU);
        vecB = bsxfun(@minus,vecB,mean(vecB,2));
        
        
        
       
        
        %vecB = imfilter(vecB,fspecial('gaussian',[1 101],31),'replicate');
        [states,prob] = hmm.Viterbi(vecB,[1 1 1]);
        
        
        
        
        % draw vertical cone lines
        coneMask = repmat(states == 3,[size(I,1) 1]);
        nonBKMask = repmat(states~=1,[size(I,1) 1]);
        out = flattenMaskOverlay(I,coneMask);
        out = I;
        imshow(out,[]);
        hold on
        
        
        
        RE2 = regionprops(logical(coneMask),'Centroid');
        RE3 = regionprops(logical(nonBKMask),'PixelIdxList');
        RE = regionprops(states~=1,'BoundingBox','PixelIdxList');
        REC = [];;
        GREEN_BOX = {};
        for r = 1:numel(RE)
            WING_MASK = zeros(1,size(I,2));
            WING_MASK(RE(r).PixelIdxList) = 1;
            WING_MASK = imdilate(WING_MASK,strel('disk',200,0));
            WING_MASK = repmat(WING_MASK,[size(I,1) 1]);
            
            
            PLANT_MASK = zeros(size(I,1),size(I,2));
            PLANT_MASK(RE3(r).PixelIdxList) = 1;
            
            RE(r).BoundingBox(2) = sidx-CROP_BUMP;
            RE(r).BoundingBox(4) = sidx2 - sidx + CROP_BUMP;
            GREEN_BOX{r} = RE(r).BoundingBox;
            rectangle('Position',RE(r).BoundingBox,'EdgeColor','g');
            tmpTOP = imcrop(I,RE(r).BoundingBox);
            
            %%% findlines
            [xy] = findHorizontalLines(tmpTOP);
            
            xy(:,1,:) = xy(:,1,:) + RE(r).BoundingBox(2);
            xy(:,2,:) = xy(:,2,:) + RE(r).BoundingBox(1);
            
            lCL = {'y' 'y'};
            for k = 1:size(xy,3)
                plot(xy(:,2,k),xy(:,1,k),'Color',lCL{k});
            end
            
            
            
            
            coneTOPY = mean(xy(:,1,k));
            PLANT_MASK(coneTOPY:end,:) = 0;
            WING_MASK(coneTOPY:end,:) = 0;
            PLANT_MASK(1:(box(2)+box(4)),:) = 0;
            WING_MASK(1:(box(2)+box(4)),:) = 0;
            
            out = flattenMaskOverlay(out,logical(PLANT_MASK),.25,'b');
            out = flattenMaskOverlay(out,logical(WING_MASK - PLANT_MASK),.05,'b');
            %PLANT_MASKWING = imdilate(PLANT_MASK,strel('line',31,0));
            %PLANT_MASKWING = PLANT_MASKWING - PLANT_MASK;
            %out = flattenMaskOverlay(out,logical(PLANT_MASKWING),.25,'b');
            
            
            REC{r} = regionprops(logical(WING_MASK),'BoundingBox');
            
            plant_Image = imcrop(I,REC{r}.BoundingBox);
            
            
            
            plant_Image_smooth = imfilter(plant_Image,fspecial('gaussian',[21],5),'replicate');
            Lab = rgb2lab(plant_Image_smooth);
            YB = (Lab(:,:,3) + 100)/200;


            GMModel = fitgmdist(YB(:),2);
            preMASK = reshape(GMModel.cluster(YB(:)),size(YB));


            [~,plantClassID] = min(GMModel.ComponentProportion);
            MASK_PLANT = preMASK == plantClassID;
            MASK_PLANT = bwareaopen(MASK_PLANT,2000);
            mi = [];
            [mi(:,1),mi(:,2)] = find(MASK_PLANT);
            mi(:,1) = mi(:,1) + REC{r}.BoundingBox(2);
            mi(:,2) = mi(:,2) + REC{r}.BoundingBox(1);
            mi = floor(mi);
            miIDX = sub2ind(size(WING_MASK),mi(:,1),mi(:,2));
            MASK_PLANT = zeros(size(WING_MASK));
            MASK_PLANT(miIDX) =  1;
            
            
            MASTER_PLANT_STACK = cat(4,MASTER_PLANT_STACK,imresize(plant_Image,[1600 1600]));
            
            TOP_STOCK = cat(4,TOP_STOCK,imresize(tmpTOP,[100 500]));
            
            out = flattenMaskOverlay(out,logical(MASK_PLANT),.6,'g');
        end
        
        
        imshow(out,[]);
        hold on
        for r = 1:numel(REC)
            rectangle('Position',REC{r}.BoundingBox,'EdgeColor','k')
        end
        
        for r = 1:numel(GREEN_BOX)
            rectangle('Position',GREEN_BOX{r},'EdgeColor','g');
        end
        
        
        %kidx = GMModel.cluster([vecB']);
        %plot(100*kidx +  3*size(I,1)/4,'k');
        CL = {'r' 'g' 'b'};
        for k = 1:size(vecB,1)
            plot(200*vecB(k,:)' + size(I,1)/4,CL{k});
        end
        
        
        plot(100*(states) +  size(I,1)/2,'m');
        
        
        hold off
        
        drawnow
        
       
        
        HMC{e} = [states;vecB];
        
        %%%% gather for conetainers
        
        
        
        
        
        
        
        if gatherQR
            %%%% start QR data gather
            Lab = rgb2lab(I);
            RED = Lab(:,:,2) > 25.5;
            RED = imclose(RED,strel('square',51));
            RED = bwlarge(RED);
            R = regionprops(RED);
            box = R(1).BoundingBox;
            box(1:2) = box(1:2) - pad;
            box(3:4) = box(3:4) + 2*pad;
            frame = imcrop(RED,box);
            frame = imclose(frame,strel('square',11));
            holes = imfill(frame,'holes') - frame;
            smallHoles = holes - bwareaopen(holes,1000);
            frame = logical(frame + smallHoles);
            qr = imcrop(I,box);






            if preD
                tmpI = imresize(qr,1/FR);
                X = im2colF(double(tmpI),[WZ 3],[1 1 1]);
                tmpX = reshape(X,[WZ 3 size(X,2)]);
                ttr = [];
                ttc = [];
                for s = 1:8
                    preY = redCornerNet{s}.predict(tmpX,'MiniBatchSize',1000);
                    pY = col2im(preY(:,2),[WZ],[size(tmpI,1) size(tmpI,2)]);
                    [~,midx] = max(pY(:));
                    [ttr(s),ttc(s)] = ind2sub(size(pY),midx);
                end
                %{
                imshow(tmpI,[]);
                hold on
                plot(ttc+10,ttr+10,'g*');
                 drawnow
                hold off
                %}

                MASTER{e} = [ttc' tte'];
                MASTERI{e} = tmpI;

                toc



            else


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                BLUE = Lab(:,:,3) < -8;
                BLUE = imclose(BLUE,strel('square',15));
                BLUE = bwlarge(BLUE,24);
                BLUEframe = imcrop(BLUE,box);
                BLUEframe = imdilate(BLUEframe,strel('square',9));
                %BLUEframe = imdilate(BLUEframe,strel('line',13,0));
                %BLUEframe = imdilate(BLUEframe,strel('line',13,90));
                %BLUEframe = imerode(BLUEframe,strel('line',13,0));
                %BLUEframe = imerode(BLUEframe,strel('line',13,90));
                %BLUEframe = imerode(BLUEframe,strel('square',5));
                holes = imfill(BLUEframe,'holes') - BLUEframe;
                smallHoles = holes - bwareaopen(holes,1000);
                BLUEframe = logical(BLUEframe + smallHoles);
                BLUEskeleton = bwmorph(BLUEframe,'skel',inf');
                BLUEskeleton = bwmorph(BLUEskeleton,'spur',inf);
                BLUEskeleton = bwmorph(BLUEskeleton,'skel',inf');
                BLUEskeleton = imdilate(BLUEskeleton,strel('square',1));
                BLUEsquares = imfill(BLUEskeleton,'holes');
                BLUEsquares = imdilate(BLUEsquares,strel('square',4));
                dB = bwboundaries(BLUEsquares);
                boxCorners = [];
                for b = 1:numel(dB)
                    K = cwtK_closed_peaks(dB{b}(1:(end-1),:),{11},50,2*10^-8);
                    boxCorners(:,:,b) = flipdim(dB{b}(K.pIDX,:),2);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                skeleton = bwmorph(frame,'skel',inf');
                skeleton = bwmorph(skeleton,'spur',inf);
                skeleton = bwmorph(skeleton,'skel',inf');
                branchPoints = bwmorph(skeleton,'branchpoints',1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% find the point orders
                br = [];
                [br(:,2),br(:,1)] = find(branchPoints==1);
                [~,TOPPOINT_IDX] = min(br(:,2));
                [~,LEFTPOINT_IDX] = min(br(:,1));
                [~,RIGHTPOINT_IDX] = max(br(:,1));
                MIDDLEPOINT_IDX = setdiff(1:4,[TOPPOINT_IDX LEFTPOINT_IDX RIGHTPOINT_IDX]);
                redPOINTS(2,:) = br(TOPPOINT_IDX,:);
                redPOINTS(4,:) = br(LEFTPOINT_IDX,:);
                redPOINTS(5,:) = br(MIDDLEPOINT_IDX,:);
                redPOINTS(6,:) = br(RIGHTPOINT_IDX,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                skeletonFrame = skeleton;
                skeletonFrame((redPOINTS(5,2) - 50):(redPOINTS(5,2) + 50),(redPOINTS(5,1) - 50):(redPOINTS(5,1) + 50)) = 0;
                skeletonFrame_lessTOP = skeletonFrame;
                skeletonFrame_lessTOP((redPOINTS(2,2) - 50):(redPOINTS(2,2) + 50),(redPOINTS(2,1) - 50):(redPOINTS(2,1) + 50)) = 0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace edge from left to top
                sk = [];
                [sk(:,2),sk(:,1)] = find(skeletonFrame);
                sourcePoint = snapTo(sk,redPOINTS(4,:));
                targetPoint = snapTo(sk,redPOINTS(2,:));
                ADJ = Radjacency(sk',2);
                [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
                path = sk(pathIDX,:);
                K = cwtK_filter(path,{7});
                K.K(1:50) = 0;
                K.K(end-50:end) = 0;
                [~,ptIDX] = max(abs(K.K));
                redPOINTS(1,:) = path(ptIDX,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace edge from top to right
                sk = [];
                [sk(:,2),sk(:,1)] = find(skeletonFrame);
                sourcePoint = snapTo(sk,redPOINTS(2,:));
                targetPoint = snapTo(sk,redPOINTS(6,:));
                ADJ = Radjacency(sk',2);
                [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
                path = sk(pathIDX,:);
                K = cwtK_filter(path,{7});
                K.K(1:50) = 0;
                K.K(end-50:end) = 0;
                [~,ptIDX] = min(K.K);
                redPOINTS(3,:) = path(ptIDX,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % trace edge from left to right
                sk = [];
                [sk(:,2),sk(:,1)] = find(skeletonFrame_lessTOP);
                sourcePoint = snapTo(sk,redPOINTS(4,:));
                targetPoint = snapTo(sk,redPOINTS(6,:));
                ADJ = Radjacency(sk',2);
                [pathIDX] = dijkstra(ADJ,sourcePoint,targetPoint);
                path = sk(pathIDX,:);
                K = cwtK_filter(path,{7});
                K.K(1:50) = 0;
                K.K(end-50:end) = 0;
                [~,ptIDX] = min(K.K);
                f1 = imdilate(K.K,ones(100,1)) == K.K;
                f2 = K.K > .02;
                fidx = find(f1.*f2);
                redPOINTS(7:8,:) = path(fidx,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




                %imshow(skeleton,[]);




                %{
                tmpQ = im2colF(qr,[21 21 3],[1 1 1]);

                if exist(qE)
                    qqC = PCA_REPROJ_T(tmpQ,qE,qU);
                    for p = 1:size(qqC,1)
                        P(:,:,p) = col2im(qqC(p,:),[21 21],[size(qr,1) size(qr,2)]);
                    end
                end

                tmpF = im2colF(double(frame),[21 21],[1 1]);
                tmpS = im2colF(double(skeleton),[21 21],[1 1]);
                kidx = find(tmpF((end-1)/2,:) == 1);
                kidxS = find(tmpS((end-1)/2,:) == 1);
                tmpQ = tmpQ(:,kidx);
                QR_PIECES = [QR_PIECES tmpQ];


                QRSHEET = [QRSHEET;imresize(qr,round([500 700]/4))];
                %}
                figure;
                imshow(qr,[]);
                hold on
                pause(.3)
                plot(br(:,1),br(:,2),'bo');
                plot(redPOINTS(:,1),redPOINTS(:,2),'g*');
                CL = {'r*','m*','k*','c*'};
                for b = 1:size(boxCorners,3)
                    for p = 1:size(boxCorners,1)
                        plot(boxCorners(p,1,b),boxCorners(p,2,b),CL{p})
                    end
                end
                blue{cnt} = boxCorners;
                clickPoints{cnt} = redPOINTS;
                %clickedFileList{cnt} = imread(FileList{e});
                QR{cnt} = qr;
                cnt = cnt + 1;

                hold off;
                drawnow

            end
        end
    catch ME
        getReport(ME)
    end
end
%% view plant images
close all
for e = 1:size(MASTER_PLANT_STACK,4)
    imshow(MASTER_PLANT_STACK(:,:,:,e),[]);
    drawnow
end
%% view plant images
close all
for e = 1:size(MASTER_PLANT_STACK,4)
    iF = imfilter(MASTER_PLANT_STACK(:,:,:,e),fspecial('gaussian',21,5),'replicate');
    
   
    Lab = rgb2lab(iF);
    YB = (Lab(:,:,3) + 100)/200;
    
    
    GMModel = fitgmdist(YB(:),2);
    preMASK = reshape(GMModel.cluster(YB(:)),size(YB));
    
    
    [~,plantClassID] = min(GMModel.ComponentProportion);
    MASK = preMASK == plantClassID;
    MASK = bwareaopen(MASK,2000);
    %tMASK = graythresh(YB);
    %MASK = YB > tMASK;
    SKEL = bwmorph(MASK,'skel',inf);
    bp = [];
    ep = [];
    [bp(:,2),bp(:,1)] = find(bwmorph(SKEL,'branchpoints'));
    [ep(:,2),ep(:,1)] = find(bwmorph(SKEL,'endpoints'));
    E = edge(Lab(:,:,2));
   
    E = imdilate(E,strel('disk',5));
    %MASK = imfill(E,'holes');
    out = flattenMaskOverlay(MASTER_PLANT_STACK(:,:,:,e),MASK,.5,'r');
    out = flattenMaskOverlay(out,SKEL,.5,'b');
    imshow(out,[]);
    hold on
    plot(bp(:,1),bp(:,2),'k*');
    plot(ep(:,1),ep(:,2),'ko');
    hold off
    drawnow
end
%% finding conetainer -- step 1 after gather start PCA for bottom container
%bottomSTACK = rgb2gray(bottomSTACK);
[bS bC bU bE bL bERR bLAM] = PCA_FIT_FULL_T(bottomSTACK,3);
GMModel = fitgmdist([bC'],3);
%% finding conetainer - step 2 make hmm 
kidx = GMModel.cluster(bC');
backGroundIDX = find(kidx==2);
coneIDX = find(kidx==3);
containerIDX = find(kidx==1);

background_node = hmm_node('background');
cone_node = hmm_node('cone');
container_node_PRE = hmm_node('container_pre');
container_node_POST = hmm_node('container_post');

backgroundCov = cov(bC(:,backGroundIDX)');
backgroundMean = mean(bC(:,backGroundIDX)');

coneCov = cov(bC(:,coneIDX)');
coneMean = mean(bC(:,coneIDX)');

conetainerCov = cov(bC(:,containerIDX)');
conetainerMean = mean(bC(:,containerIDX)');


backgroundDis = myProb(backgroundMean,backgroundCov);
coneDis = myProb(coneMean,coneCov);
conetainerDis = myProb(conetainerMean,conetainerCov);
    
background_node.attachDistribution(backgroundDis,1);
cone_node.attachDistribution(coneDis,1);
container_node_PRE.attachDistribution(conetainerDis,1);
container_node_POST.attachDistribution(conetainerDis,1);

backGroundTobackGround = constantTransitionFunction(.5);
backGroundTobackGround = VarheavisideTransitionFunction(200,@(x,y)lt(x,y),.9,1);
backGroundTocontainer = constantTransitionFunction(.001);
backGroundTocontainer = VarheavisideTransitionFunction(100,@(x,y)ge(x,y),.9,1);

containerTocontainer = constantTransitionFunction(.8);
containerTocontainer = VarheavisideTransitionFunction(400,@(x,y)lt(x,y),.9,1);
%containerTocontainer = heavisideTransitionFunction(550,@(x,y)lt(x,y));
%containerTocontainer = heavisideTransitionFunction(300,@(x,y)lt(x,y));
%containerPRETocontainerPOST = constantTransitionFunction(.0001);
containerTobackground = constantTransitionFunction(.1);
%containerTobackground = VarheavisideTransitionFunction(400,@(x,y)lt(x,y),.9999,1);
%containerTobackground = heavisideTransitionFunction(550,@(x,y)ge(x,y));
%containerTobackground = heavisideTransitionFunction(300,@(x,y)ge(x,y));
%containerTocone = constantTransitionFunction(.00000001);
containerTocone = VarheavisideTransitionFunction(400,@(x,y)lt(x,y),.01,1);
%containerTocone = heavisideTransitionFunction(550,@(x,y)ge(x,y));

coneTocone = constantTransitionFunction(.9);
coneTocone = heavisideTransitionFunction(215,@(x,y)lt(x,y));
coneTocone = VarheavisideTransitionFunction(215,@(x,y)lt(x,y),1,1);
coneToconetainer = constantTransitionFunction(.001);
coneToconetainer = VarheavisideTransitionFunction(215,@(x,y)ge(x,y),1,1);
%coneToconetainer = heavisideTransitionFunction(150,@(x,y)ge(x,y));

background_node.attachNode(background_node,backGroundTobackGround);
background_node.attachNode(container_node_PRE,backGroundTocontainer);

container_node_PRE.attachNode(container_node_PRE,containerTocontainer);
container_node_PRE.attachNode(cone_node,containerTocone);
%container_node_PRE.attachNode(container_node_POST,containerPRETocontainerPOST);

container_node_POST.attachNode(container_node_POST,containerTocontainer);
container_node_POST.attachNode(background_node,containerTobackground);

cone_node.attachNode(cone_node,coneTocone);
cone_node.attachNode(container_node_POST,coneToconetainer);


hmm = my_hmm();
hmm.addNode(background_node);
hmm.addNode(container_node_PRE);
hmm.addNode(cone_node);
hmm.addNode(container_node_POST);

hmm.dn = ones(size(I,2),1);
%%
hmm.update(HMC,[1 1 1]);
%% gather for PC analysis
close all
gTMP = [];
gTMP2 = [];
gTMP3 = [];
for e = 1:size(TOP_STOCK,4)
    imshow(TOP_STOCK(:,:,:,e),[]);
    title(e)
    hold on
    
    %pause(.3)
    gr = rgb2gray(TOP_STOCK(:,:,:,e));
    [d1 d2] = gradient(gr);
    BW = edge(gr);
    
    findHorizontalLines(I)
    
    [H,T,R] = hough(BW','Theta',linspace(-10,10,100));
    P  = houghpeaks(H,1);
    lines = houghlines(BW',T,R,P,'FillGap',200,'MinLength',50);
    
    for k = 1:numel(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,2),xy(:,1),'LineWidth',2,'Color','green');
        
    end
    
    CHOP = mean(xy(:,1));
    mBW = BW(1:(CHOP-2),:);
    
    [H,T,R] = hough(mBW','Theta',linspace(-10,10,100));
    P  = houghpeaks(H,1);
    lines = houghlines(mBW',T,R,P,'FillGap',200,'MinLength',50);
    
    for k = 1:numel(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,2),xy(:,1),'LineWidth',2,'Color','red');
        
    end
    
    
    
    
    gTMP3 = [gTMP3 mean(gr,2)];
    
    gTMP = [gTMP bsxfun(@minus,gr,mean(d2,2))];
    gTMP2 = [gTMP2;bsxfun(@minus,gr,mean(d2,1))];
    drawnow
    hold off
    e
    waitforbuttonpress
    
end
%% perform PC analysis
[cS cC cU cE cL cERR cLAM] = PCA_FIT_FULL_T(gTMP,3);
[cS2 cC2 cU2 cE2 cL2 cERR2 cLAM2] = PCA_FIT_FULL(gTMP2,3);
[cS3 cC3 cU3 cE3 cL3 cERR3 cLAM3] = PCA_FIT_FULL_T(gTMP3,3);
GMModel_ConeLine = fitgmdist([gTMP3(:)],3);

%% view PC analysis
close all
gTMP = [];
gTMP2 = [];
for e = 1:size(TOP_STOCK,4)
    gr = rgb2gray(TOP_STOCK(:,:,:,e));
    tmpC1 = PCA_REPROJ_T(bsxfun(@minus,gr,mean(gr,1)),cE,cU);
    tmpC2 = PCA_REPROJ(bsxfun(@minus,gr,mean(gr,2)),cE2,cU2);
    
    tmpC3 = [mean(gr,2)];
    kline = GMModel_ConeLine.cluster(tmpC3);
    lidx = find(kline == 2);
    lidx = lidx(end);
    
    imshow(TOP_STOCK(:,:,:,e),[]);
    hold on
    plot(1:size(TOP_STOCK,2),lidx*ones(size(TOP_STOCK,2)),'r');
    hold off
    drawnow
    %{
    hold on
    for k = 1:size(tmpC1,1)
        plot(tmpC1(k,:)'*20 + size(TOP_STOCK,1)/2,CL{k});
    end
    plot(40*tmpC2'+size(TOP_STOCK,2)/2,1:size(TOP_STOCK,1));
    drawnow
    waitforbuttonpress
    %}
end
%% gather QR for croppping
autoPath_good = '/mnt/tetra/nate/autoRCNN/good/';
autoPath_bad = '/mnt/tetra/nate/autoRCNN/bad/';
close all
pad = 50;
QRSHEET = [];
QR_PIECES = [];

cnt = 1;
clickedFileList = {};
clickPoints = {};
preD = false;
trainTable = table;
bottomSTACK = [];
conetainerBUMP = 150;
for e = 1:100%:1:7000%numel(FileList)
    
        [pth,nm,ext] = fileparts(FileList{e});
    
        tic
        I = double(imread(FileList{e}))/255;
        I(:,(end-100):end,:) = [];

        
        
        

        Lab = rgb2lab(I);
        RED = Lab(:,:,2) > 25.5;
        RED = imclose(RED,strel('square',51));
        RED = bwlarge(RED);
        R = regionprops(RED);
        box = R(1).BoundingBox;
        box(1:2) = box(1:2) - pad;
        box(3:4) = box(3:4) + 2*pad;
        frame = imcrop(RED,box);
        frame = imclose(frame,strel('square',11));
        holes = imfill(frame,'holes') - frame;
        smallHoles = holes - bwareaopen(holes,1000);
        frame = logical(frame + smallHoles);
        qr = imcrop(I,box);
        
        
        thumb = imresize(I,.15);
        
        
        EN = entropyfilt(rgb2gray(thumb));
        BK = EN < 2.3;
        
        BK_thumb = bsxfun(@times,BK,thumb);
        fidx = find(~BK);
        for k = 1:3
            tmp = BK_thumb(:,:,k);
            tmp(fidx) = NaN;
            BK_thumb(:,:,k) = tmp;
        end
        BK_thumb1 = nanmean(BK_thumb,1);
        BK_thumb2 = nanmean(BK_thumb,2);
        BK_thumb1 = repmat(BK_thumb1,[size(BK_thumb,1) 1]);
        BK_thumb2 = repmat(BK_thumb2,[1 size(BK_thumb,2)]);
        BK_thumb1 = imfilter(BK_thumb1,fspecial('disk',31),'replicate');
        BK_thumb2 = imfilter(BK_thumb2,fspecial('disk',31),'replicate');
        BK_thumbT = .5*(BK_thumb1 + BK_thumb2);
        
        
        
        
        IGONE = thumb;
        box = round(box*.15);
        IGONE(box(2):(box(2)+box(4)),box(1):(box(1)+box(3)),:) = BK_thumbT(box(2):(box(2)+box(4)),box(1):(box(1)+box(3)),:);
        
        QR = insertObjectAnnotation(thumb, 'rectangle', box, 'QR');
        imshow(QR,[]);
        drawnow
        
        trainTable(e,'fileName') = {[autoPath_good nm '.tif']};
        trainTable(e,'QR_sheet') = {{[box]}};
        imwrite(thumb,[autoPath_good nm '.tif']);
        imwrite(IGONE,[autoPath_bad nm '.tif']);
        
        
    
end
%%
negativeImages = imageDatastore(autoPath_bad);
trainCascadeObjectDetector('QR.xml',trainTable, ...
    negativeFolder,'FalseAlarmRate',0.02,'NumCascadeStages',5,'FeatureType','HOG','NegativeSamplesFactor',50);
%%
close all
detector = vision.CascadeObjectDetector('QR.xml');
bx = step(detector,thumb);
IFaces = insertObjectAnnotation(thumb, 'rectangle', bx, 'Face');
imshow(IFaces)
%%
CROPlayers = [imageInputLayer([28 28 3])
        convolution2dLayer([5 5],10)
        reluLayer()
        maxPooling2dLayer(3,'Stride',2);
        fullyConnectedLayer(2)
        softmaxLayer()
        classificationLayer()];
    

options = trainingOptions('sgdm','Plots','training-progress','MaxEpochs',10,'InitialLearnRate',0.0001,'ExecutionEnvironment','auto');

    detector = trainRCNNObjectDetector(trainTable,CROPlayers,options);
%%
close all
pad = 50;
QRSHEET = [];
QR_PIECES = [];

cnt = 1;
clickedFileList = {};
clickPoints = {};
preD = true;

parfor e = 1:1:7000%numel(FileList)
    try
        tic
        I = double(imread(FileList{e}))/255;
        I(:,(end-100):end,:) = [];


        Lab = rgb2lab(I);
        RED = Lab(:,:,2) > 25.5;
        RED = imclose(RED,strel('square',51));
        RED = bwlarge(RED);
        R = regionprops(RED);
        box = R(1).BoundingBox;
        box(1:2) = box(1:2) - pad;
        box(3:4) = box(3:4) + 2*pad;
        frame = imcrop(RED,box);
        frame = imclose(frame,strel('square',11));
        holes = imfill(frame,'holes') - frame;
        smallHoles = holes - bwareaopen(holes,1000);
        frame = logical(frame + smallHoles);
        qr = imcrop(I,box);

        
        
    
        if preD
            tmpI = imresize(qr,1/FR);
            X = im2colF(double(tmpI),[WZ 3],[1 1 1]);
            tmpX = reshape(X,[WZ 3 size(X,2)]);
            ttr = [];
            ttc = [];
            for s = 1:8
                preY = redCornerNet{s}.predict(tmpX,'MiniBatchSize',1000);
                pY = col2im(preY(:,2),[WZ],[size(tmpI,1) size(tmpI,2)]);
                [~,midx] = max(pY(:));
                [ttr(s),ttc(s)] = ind2sub(size(pY),midx);
            end
            %{
            imshow(tmpI,[]);
            hold on
            plot(ttc+10,ttr+10,'g*');
             drawnow
            hold off
            %}
            
            MASTER{e} = [ttc' ttr'];
            MASTERI{e} = tmpI;
           
            toc
        else
        end
    catch ME
        getReport(ME)
        
    end
end
%%
%[X,map] = rgb2ind(QR{1},3,'dither');
WZ = [21 21];
YS = [];
XS = [];
FR = 4;

M = 100;
% select one byte of information at random
pixel_select = randi(prod([WZ 3]),8,M);
bit_select = randi(8,8,M);
MX = [];
preD = false;
for e = 1:numel(QR)
    tmpI = QR{e};
    tmpI = imresize(tmpI,1/FR);
    
    X = im2colF(double(tmpI),[WZ 3],[1 1 1]);
    Y = [];
    
    tmpX = reshape(X,[WZ 3 size(X,2)]);
    
    
    
    if preD
        for s = 1:8
            preY = redCornerNet{s}.predict(tmpX,'MiniBatchSize',1000);
            pY = col2im(preY(:,2),[WZ],[size(tmpI,1) size(tmpI,2)]);
            [~,midx] = max(pY(:));
            [ttr(s),ttc(s)] = ind2sub(size(pY),midx);
        end
        imshow(tmpI,[]);
        hold on
        plot(ttc+10,ttr+10,'g*');
    end
    
    parfor pt = 1:8
        tmp = zeros(size(tmpI,1),size(tmpI,2));
        tmp(round(clickPoints{e}(pt,2)/FR),round(clickPoints{e}(pt,1)/FR)) = 1;
        tmp = imdilate(tmp,strel('disk',2,0));
        tmpY = im2colF(tmp,[WZ],[1 1]);
        Y(pt,:) = pt*tmpY((end-1)/2,:);
        %out = flattenMaskOverlay((tmpI+.1)/1.5,logical(tmp));
        %imshow(out,[]);
        %drawnow
        %pt
    end
    Y = sum(Y,1);
    %{
    Lab = rgb2lab(tmpI);
    tmpI = bindVec(Lab(:,:,2));
    tmpI = tmpI * 255;
    %}
    
   
    %J = autoenc.encode(double(X));
    %MX = [MX,double(X)];
    
  
    
    
   
    

    fidx1 = find(Y~=0);
    fidx0 = find(Y==0);


    fidx0 = fidx0(randi(numel(fidx0),1,100));
    fidxT = [fidx0';fidx1'];
    
    YS = [YS Y(fidxT)];
    XS = [XS X(:,fidxT)];
    
    e
end
%%
X = reshape(XS,[WZ 3 size(XS,2)]);
layers = [imageInputLayer([size(X,1) size(X,2) size(X,3)],'Normalization','none');
          convolution2dLayer(5,5);
          reluLayer();
          maxPooling2dLayer(3,'Stride',2);
          convolution2dLayer(5,5);
          reluLayer();
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];

      
      
%{
X = X(:,:,:,cidx);
YS = YS(cidx);
%}
for e = 1:8
    fidx1 = find(YS==e);
    fidx0 = find(YS~=e);
    
    holdOutNumber1 = round(.80*numel(fidx1));
    holdOutNumber0 = round(.80*numel(fidx0));
    
    fidx1 = fidx1(randperm(numel(fidx1)));
    fidx0 = fidx0(randperm(numel(fidx0)));
    
    TRAIN1 = fidx1(1:holdOutNumber1);
    TEST1 = fidx1(holdOutNumber1:end);
    TRAIN0 = fidx0(1:holdOutNumber0);
    TEST0 = fidx0(holdOutNumber0:end);
    
    trainX = X(:,:,:,[TRAIN1 TRAIN0]);
    trainY = YS([TRAIN1 TRAIN0]);
    
    testX = X(:,:,:,[TEST1 TEST0]);
    testY = YS([TEST1 TEST0]);
   

    options = trainingOptions('sgdm','Shuffle','every-epoch','ValidationFrequency',300,'ValidationData',{testX,categorical(double(testY'==e))},'Plots','training-progress','MaxEpochs',10,'InitialLearnRate',0.01,'ExecutionEnvironment','parallel');

    redCornerNet{e} = trainNetwork(trainX,categorical(double(trainY'==e)),layers,options);
end
%%

%%
hiddenSize = 10;
autoenc = trainAutoencoder(MX,hiddenSize,...
        'EncoderTransferFunction','logsig',...
        'DecoderTransferFunction','purelin',...
        'L2WeightRegularization',0.01,...
        'SparsityRegularization',4,...
        'SparsityProportion',0.03);

%%
for r = 1:size(XS,2)

    close all
    iterations = 200;
    repeats = 1;
    MAX_B = 9;
    
    BYTE_SIZE = 3;
    MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide

    MAX_P_NUMBER = min(255,2^MAX_B);
    qXS = quantumEncode(MAX_B,XS(:,r),true);
    evalProgram = @(X)evalSimMetrics(X,uint8(YS),'hamming');
    [STORE,DIS] = findProgram(1,[figure;figure],repeats,iterations,evalProgram,qXS,YS,100,MAX_WID,BYTE_SIZE,MAX_P_NUMBER);

    [~,pidx] = min(DIS(:,1));
    TR = boo(STORE(pidx,:),qXS);
    
    
    DIS_STORE{r} = DIS;
    PROG_STORE{r} = STORE;
end
%%
close all
for e = 1
    tmpI = QR{e};
    tmpI = imresize(tmpI,1/FR);
    tmp = zeros(size(tmpI,1),size(tmpI,2));
    tmp(round(clickPoints{e}(1,2)/FR),round(clickPoints{e}(1,1)/FR)) = 1;
    tmp = imdilate(tmp,strel('disk',5,0));
    X = rgb2ind(tmpI, map,'dither');
    
    
    
    Lab = rgb2lab(tmpI);
    tmpI = bindVec(Lab(:,:,2));
    X = tmpI * 255;
    
    
    
    
    X = uint8(im2colF(double(X),WZ,[1 1]));


    %[X] = squeezeQbits(X(sidx,:)',2,4,'uint16');
    
    
    [X] = squeezeQbits2(X(sidx,:)',{[8],[8],[8],[8],[8],[8],[8],[8],[8]},4,'uint16');
    
    
    
    X = quantumEncode(MAX_B,X);
    preY = boo(STORE(pidx,:),X);
    preY = col2im(preY,WZ,size(tmpI),'sliding');
    %preY = imdilate(preY,strel('disk',2,0));
    
    R = regionprops(preY==1,'Centroid');
    figure;
    imshow(preY==1,[]);
    
    figure;
    imshow(tmpI,[]);
    pause(1)
    hold on
    for p = 1:numel(R)
        plot(R(p).Centroid(1)+10,R(p).Centroid(2)+10,'r*');
    end
    drawnow
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make noise images for red sparkle - ver2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
SPboxSequence = [1700 1700];
SPzoomSequence = [.15];

cnt = 1;
NEWY = [];
NEWX = [];
MAHGY = [];
CENTER = [];
for e = 1:100%numel(QR)
    for m = 1:1
        mag = (.8*rand(1)+.6);
        mag = 1;
        
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = double(QR{e});
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        nIMGx = rand(size(nIMGx));
        
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k));
        end
        
        
        GOOD_red_corners = clickPoints{e}*SPzoomSequence_tmp;
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        
        
        
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        drawnow
        hold off
        
        
        
        
        
        
        NEWX(:,:,:,cnt) = nIMGx;
        MAHGY(cnt) = mag;
        CENTER(:,cnt) = mean(GOOD_red_corners,1);
        NEWY(:,:,cnt) = GOOD_red_corners;
        
        
        
        
        
        
        %GOOD_red_corners = bsxfun(@minus,GOOD_red_corners,mean(GOOD_red_corners,1));
        %NEWY(cnt,:) = GOOD_red_corners(:);%*SPzoomSequence^-1;
        cnt = cnt + 1;
    end
    e
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make noise images for red sparkle - ver2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
SPboxSequence = [1300 1300];
SPzoomSequence = [.15];
cnt = 1;
NEWY = [];
NEWX = [];
window_size = [15 15 3];
dilateAmount = 3;
sampleN_0 = 200;



for e = 1:100%numel(QR)
    
    for m = 1:5
        mag = (.8*rand(1)+.6);
        mag = 1;
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = double(QR{e});
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        
        
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k));
        end
        
        
        GOOD_red_corners = clickPoints{e}*SPzoomSequence_tmp;
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        
        
        
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        drawnow
        hold off
        
        
        tmpF = im2colF(nIMGx,window_size,[1 1 1]);
        
        class_label = [];
        
        for class = 1:size(GOOD_red_corners,1)
            tmp_mask = zeros(size(nIMGx,1),size(nIMGx,2));
            tmp_mask(round(GOOD_red_corners(class,2)),round(GOOD_red_corners(class,1))) = 1;
            tmp_mask = imdilate(tmp_mask,strel('disk',dilateAmount,0));
            tmp_mask = im2colF(tmp_mask,window_size(1:2),[1 1]);
            tmp_mask = tmp_mask((end-1)/2,:);
            class_label(class,:) = tmp_mask;
            
            
            fidx1 = find(tmp_mask==1);
            fidx0 = find(tmp_mask==0);
            fidx0 = fidx0(randperm(numel(fidx0)));
            
            
            samp1 = tmpF(:,fidx1);
            samp0 = tmpF(:,fidx0);
            
            
            if e == 1
                Xdata{class} = [];
                Ydata{class} = [];
            end
            
            Xdata{class} = [Xdata{class},samp1,samp0(:,1:sampleN_0)];
            Ydata{class} = [Ydata{class},ones(1,numel(fidx1)),zeros(1,sampleN_0)];
            
            
        end
        
        
        %NEWX(:,:,:,cnt) = nIMGx;
        %NEWY(:,:,cnt) = GOOD_red_corners;
        %GOOD_red_corners = bsxfun(@minus,GOOD_red_corners,mean(GOOD_red_corners,1));
        %NEWY(cnt,:) = GOOD_red_corners(:);%*SPzoomSequence^-1;
        cnt = cnt + 1;
    end
    e
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for class = 1:numel(Xdata)
    sz = size(Xdata{class});
    Xdata{class} = reshape(Xdata{class},[window_size sz(end)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% class approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
net = {};
for class = 1:numel(Xdata)
    close all
    options = trainingOptions('sgdm',...
        'InitialLearnRate',.01,...
        'MaxEpochs',10,...
        'Plots','training-progress',...
        'ExecutionEnvironment','parallel',...
        'Momentum',.9);
    layers = [
        imageInputLayer(window_size,'Normalization','None')
        convolution2dLayer(5,8,'Padding','same')
        reluLayer()
        fullyConnectedLayer(2)
        softmaxLayer()
        classificationLayer()
        ];
    net{class} = trainNetwork(Xdata{class},categorical(Ydata{class}'),layers,options);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% class approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
for e = 1000%numel(QR)
    
        mag = (.8*rand(1)+.6);
        mag = 1;
        SPzoomSequence_tmp = mag*SPzoomSequence;
       
        I = double(QR{e});
        
        I = imresize(I,SPzoomSequence_tmp);
        Isz = size(I);
        
        nIMGx = rand([fliplr(SPboxSequence) 3]);
        nIMGx = imresize(nIMGx,SPzoomSequence);
        
        
        %DIS = round(size(nIMGx)/2 - Isz/2);
        DIS(2) = size(nIMGx,2) - size(I,2);
        DIS(1) = size(nIMGx,1) - size(I,1);

        POSD = [randi(DIS(2),1) randi(DIS(1),1)];
        for k = 1:3
            nIMGx(POSD(2):(POSD(2)+size(I,1)-1),POSD(1):(POSD(1)+size(I,2)-1),k) = double(I(:,:,k));
        end
        
        GOOD_red_corners = clickPoints{e}*SPzoomSequence_tmp;
        GOOD_red_corners = bsxfun(@plus,GOOD_red_corners,POSD);
        
        tmpF = im2colF(nIMGx,window_size,[1 1 1]);
        tmpF = reshape(tmpF,[window_size size(tmpF,2)]);
        M = [];
        for class = 1:numel(net)
            preY = net{class}.predict(tmpF);
            preY = preY(:,2);
            [~,midx] = max(preY);
            preY = col2im(preY,window_size(1:2),[size(nIMGx,1) size(nIMGx,2)]);
            [M(class,2),M(class,1)] = ind2sub(size(preY),midx);
           
        end
        M = bsxfun(@plus,M,(window_size(1:2)-1)/2);
        
        imshow(nIMGx,[]);
        hold on
        plot(GOOD_red_corners(:,1),GOOD_red_corners(:,2),'r*')
        plot(M(:,1),M(:,2),'go')
        drawnow
        hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train regress center of mass and resize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
options = trainingOptions('sgdm',...
    'InitialLearnRate',.0000001,...
    'MaxEpochs',2000,...
    'Plots','training-progress',...
    'ExecutionEnvironment','parallel',...
    'Momentum',.9);
szX = size(NEWX);
layers = [
    imageInputLayer(szX(1:2),'Normalization','None')
    convolution2dLayer(21,5)
    reluLayer()
    fullyConnectedLayer(3)
    regressionLayer()];
Y = [CENTER' MAHGY'];
%[nY,uY,sY] =  zscore(Y);
MAG_LOC_net = trainNetwork(NEWX(:,:,1,:),Y,layers,options);
%MAG_LOC_net = trainNetwork(NEWX(:,:,1,:),Y,MAG_LOC_net.Layers,options);
%%
%nnY = bsxfun(@minus,Y,uY);
%nnY = bsxfun(@times,nnY,sY.^-1);
MAG_LOC_net.predict(NEWX(:,:,1,1))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate grid for scatter zoom for pots
rawI = double(imread(FileList{e}))/255;
%%
rI = imresize(rawI,.15);
%rI = imresize(rI,.7);
%%
close all

[potY,potSZ] = genIXgrid2(size(rI),[size(NEWX,1)/4 size(NEWX,2)/4],[0 0],[(size(NEWX,1)-1)/2 (size(NEWX,2)-1)/2]);

imshow(rI,[]);hold on;
plot(potY(:,1),potY(:,2),'g*')
close all
N = 5;
boxSequence = repmat(boxSequence,[N 1]);
%potY = flipdim(potY,2);
[potY] = runZoomSequence2(rI,MAG_LOC_net,potY,boxSequence,ones(1,N),uY,sY,ones(1,N));
potY = flipdim(potY,2);
close all
imshow(rI,[]);
hold on
for path = 1:size(potY,1)
    plot(squeeze(potY(path,1,:)),squeeze(potY(path,2,:)),'k')
end
plot(potY(:,1,1),potY(:,2,1),'g*')
plot(potY(:,1,end),potY(:,2,end),'r*')


