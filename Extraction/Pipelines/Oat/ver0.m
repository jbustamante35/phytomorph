FilePath = '/mnt/spaldingdata/nate/pan/';
FileList = {};
FileExt = {'jpg'};
FileList = fdig(FilePath,FileList,FileExt,1);
%%
e =1;
I = double(imread(FileList{e}))/255;
hsv = rgb2hsv(I);
lab = rgb2lab(I);
S = std(I,1,3);
th = graythresh(S);
M = S > th;
%M = bwlarge(M);
M = connectPlant2(M,300);
M = logical(M);
skel = bwmorph(M,'skeleton',inf);
R = regionprops(M,'Centroid','Area','Perimeter','ConvexArea','ConvexHull','FilledArea');
I2 = sum(M,2);
I1 = sum(M,1);
%dbm = I2*I1;
fidx2 = find(I2);
fidx1 = find(I1);
minLevelY = fidx2(1)*ones(1,size(I,2));
maxLevelY = fidx2(end)*ones(1,size(I,2));
minLevelX = fidx1(1)*ones(1,size(I,1));
maxLevelX = fidx1(end)*ones(1,size(I,1));
%%
close all
out = flattenMaskOverlay(I,M);
%out = flattenMaskOverlay(out,logical(dbm> 100),.07,'b');
imshow(out,[]);
hold on
plot(I1,'r')
plot(I2,1:numel(I2),'r');
plot(1:size(I,2),minLevelY,'g');
plot(1:size(I,2),maxLevelY,'g');
plot(minLevelX,1:size(I,1),'g');
plot(maxLevelX,1:size(I,1),'g');
plot(R.ConvexHull(:,1),R.ConvexHull(:,2),'b');
