FilePath = '/mnt/spaldingdata/Drone_Imagery/Arlington_2018/';
FileList = {};
FileExt = {'jpg','JPG'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
FileList = FileList(randperm(numel(FileList)));
%%
for e = 1:10
    I = double(imread(FileList{e}));
    NDVI = (I(:,:,1) - I(:,:,3)) ./ (I(:,:,1) + I(:,:,3));
    imshow(NDVI,[])
    drawnow
end
%%
[multiScaleDomain] = generateMultiResSampleDomain([100 100],[.1 .5 1],[30 30]);
[multiScaleDomain] = generateMultiResSampleDomain_DISK(300,[.1 .5 1],[50 100]);
%%

NS =  800;
%
I = double(imread(FileList{1}));
%NDVI = (I(:,:,1) - I(:,:,3)) ./ (I(:,:,1) + I(:,:,3));
NDVI = I(:,:,1);
px = randi(size(I,1)-200,100) + 100;
py = randi(size(I,1)-200,100) + 100;
PX = [px py];
[sampleData] = sampleApply(NDVI,PX,multiScaleDomain,[]);
sampleData = bsxfun(@minus,sampleData,mean(sampleData,1));
sampleData = abs(fft(sampleData,[],1));

MS = zeros([size(sampleData) NS]);
    %

parfor e = 1:NS
    I = double(imread(FileList{e}));
    %NDVI = (I(:,:,1) - I(:,:,3)) ./ (I(:,:,1) + I(:,:,3));
    NDVI = I(:,:,1);
    [d1,d2] = gradient(NDVI);
    NDVI = (d1.^2 + d2.^2).^.5;
    px = randi(size(I,1)-200,100) + 100;
    py = randi(size(I,1)-200,100) + 100;
    PX = [px py];
    [sampleData] = sampleApply(NDVI,PX,multiScaleDomain,[]);
    sampleData = bsxfun(@minus,sampleData,mean(sampleData,1));
    sampleData = abs(fft(sampleData,[],1));
    
    MS(:,:,:,:,e) = sampleData;
    %{
    imshow(NDVI,[])
    drawnow
    %}
    e
end
SZ = size(MS);
MS = reshape(MS,[SZ(1:3) prod(SZ(4:5))]);
%%
numC = 7;
SZ = size(MS);
MSx = reshape(MS,[prod(SZ(1:3)) SZ(4)]);
rm = any(isinf(MSx) |isnan(MSx),1);
MSx(:,rm) = [];
[U,E,L] = PCA_FIT_FULL_Tws(MSx,numC);
%%
numL = 15;
tC = PCA_REPROJ_T(MSx,E,U);
GMModel = fitgmdist(tC',numL);
%%
I = double(imread(FileList{400}));
I = double(imread(FileList{7}));
NDVI = I(:,:,1);

[~,BOX] = imcrop(NDVI/255);
NDVI = imcrop(NDVI,BOX);
%%
rawView = NDVI;
[d1,d2] = gradient(NDVI);
gNDVI = (d1.^2 + d2.^2).^.5;
%gNDVI = NDVI;
toView = gNDVI;
for r = 1:4
    rawView(1:99,:) = [];
    toView(1:99,:) = [];
    toView = imrotate(toView,90);
    rawView = imrotate(rawView,90);
end
[n1 n2] = ndgrid(100:5:(size(NDVI,1)-100+1),100:5:(size(NDVI,2)-100+1));
N = [n2(:) n1(:)];
tic
[sampleData] = sampleApply(gNDVI,N,multiScaleDomain,[]);
toc
sampleDataI = bsxfun(@minus,sampleData,mean(sampleData,1));
sampleDataI = abs(fft(sampleDataI,[],1));
SZ = size(sampleDataI);
sampleDataI = reshape(sampleDataI,[prod(SZ(1:3)) SZ(4)]);

%%

C = PCA_REPROJ_T(sampleDataI,E,U);
C = C';
C = reshape(C,[size(n1) numC]);
nC = [];
for k = 1:size(C,3)
    nC(:,:,k) = imresize(C(:,:,k),size(toView));
end
SZ = size(nC);
tC = reshape(nC,[prod(SZ(1:2)) SZ(3)]);
kidx = GMModel.cluster(tC);
kidx = reshape(kidx,SZ(1:2));
for k = 1:size(nC,3)
    tmp = nC(:,:,k);
    tmp = bindVec(tmp(:));
    nC(:,:,k) = reshape(tmp,size(toView,1),size(toView,2));
end
close all
%%
close all
figure;
imshow(rawView,[]);
figure;
imshow(label2rgb(uint8(kidx)),[])
vD = [1:3];
figure
imshow(nC(:,:,vD),[]);
figure;
imshow(toView,[]);
figure;
for e = 1:size(C,3)
    imshow(nC(:,:,e),[]);
    
    drawnow
    waitforbuttonpress
end

for L = 1:numL
    msk = uint8(kidx) == L;
    out = flattenMaskOverlay(bindVec(rawView),msk);
    imshow(out,[]);
    drawnow
    waitforbuttonpress
end

