I = imread('/mnt/spaldingdata/nate/carrot003.tif');
horEdgeTrim = 400;

I(:,1:horEdgeTrim,:) = [];
I = flip(I,2);
I(:,1:horEdgeTrim,:) = [];
I = flip(I,2);

G = rgb2gray(I);
sG = std(double(G),1,2);

M1 = sG/255 < th;
M1 = sG/255 < th;
M2 = repmat(M1,[1 size(G,2)]);
M2 = bwlarge(M2);
R = regionprops(logical(M2));
subG = imcrop(G,R(1).BoundingBox);
subI = imcrop(I,R(1).BoundingBox);
imshow(subG,[]);
figure;
imshow(I,[]);
out = flattenMaskOverlay(I,M2);
figure;
imshow(out,[]);
%%

%%
para.scales.value = 1;
para.resize.value = .25;
K = surKur(subG,para);
%%
close all
rootTH = 300;
bK = bindVec(K(:,:,2));
rootMask = bK < graythresh(bK);
rootMask = bwareaopen(rootMask,rootTH);
oMask = rootMask;
keepVec = rootMask(1,:);
for r = 1:4
    rootMask(:,1) = 0;
    rootMask = imrotate(rootMask,90);
end
rootMask(1,:) = keepVec;
toKeep = imclearborder(rootMask);
rootMask = (oMask == 1) & (toKeep == 0);
out = flattenMaskOverlay(subI,rootMask);
imshow(out,[]);
%%



