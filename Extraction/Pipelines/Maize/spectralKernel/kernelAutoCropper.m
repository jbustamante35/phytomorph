%%
tPath = '/mnt/snapper/nate/Hyperspectral/';
tFile = 'Image1_eup_2_FX17.hdr';
[wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo(tPath, tFile);
iFile = [tPath 'Image1_eup_2_FX17.raw'];
fid = fopen(iFile);
samplePart = readPart(fid, spatial, spectral, frames, 12);
%% view the spetral band(s)
for e = 1:size(samplePart,3)
    imshow(samplePart(:,:,e),[]);
    drawnow
end
%%
SZ = size(samplePart);
tmp = reshape(samplePart,[prod(SZ(1:2)) SZ(3)]);
toComp2 = 5;
[S2,C2,U2,E2,L2,ERR2,LAM2] = PCA_FIT_FULL(double(tmp),toComp2);
C2 = reshape(C2,[SZ(1:2) toComp2]);
%% view score image
imshow(C2(:,:,3),[]);
%% plot the eigen curves
plot(E2)
%% sweep
[sweepD] = sweepPCA(C2,E2,U2,std(C2,1,1),1:5,10);
for c = 1:size(sweepD,1)
    close all
    for e = 1:size(sweepD,2)
        plot(squeeze(sweepD(c,e,:)),'b')
        hold on
        waitforbuttonpress
    end
    waitforbuttonpress
end

%% segment out the background
close all
imshow(C2(:,:,1),[]);
drawnow
gMASK = bindVec(C2(:,:,1));
MASK = gMASK > graythresh(gMASK);
eMASK = imerode(MASK,strel('disk',5,0));
kMASK = imclearborder(eMASK);
imshow(kMASK,[]);
drawnow
%% segment out the grid based on third score of PCA which seems to highlight the spectral signature of the grid
close all
gMASK = bindVec(C2(:,:,3));
%imshow(gMASK,[]);
gridThresh = graythresh(gMASK);
gridThresh = .3;
gMASK = gMASK < gridThresh;
gMASK = bwlarge(gMASK);
innerBOX_MASK = ~gMASK;
imshow(gMASK,[]);


%% draw boxes around the seed size centers using the outer product and some other goodness
close all
s1 = sum(kMASK,1);
s2 = sum(kMASK,2);
s1 = imfilter(s1,fspecial('gaussian',[1 50],20));
s2 = imfilter(s2,fspecial('gaussian',[50 1],20));
numP1 = 5;
numP2 = 24;
CS2 = 100;
CS1 = 70;
peak1 = s1 == imdilate(s1,strel('disk',CS1,0));
peak2 = s2 == imdilate(s2,strel('disk',CS2,0));
plot(bindVec(s1),'b');
hold on
plot(peak1,'r')
close all
%waitforbuttonpress
plot(bindVec(s2));
hold on
plot(peak2,'r');
pidx1 = find(peak1);
pidx2 = find(peak2);
values1 = s1(pidx1);
values2 = s2(pidx2);
[~,sidx1] = sort(values1,'descend');
[~,sidx2] = sort(values2,'descend');
pidx1 = pidx1(sidx1);
pidx2 = pidx2(sidx2);
pidx1 = pidx1(1:numP1);
pidx2 = pidx2(1:numP2);
CM1 = zeros(1,size(kMASK,2));
CM2 = zeros(size(kMASK,1),1);
CM1(pidx1) = 1;
CM2(pidx2) = 1;
CM = CM2*CM1;
[cidx1,cidx2] = find(CM);
close all
imshow(bindVec(C2(:,:,1)),[]);
hold on
plot(cidx2,cidx1,'r*');
BOXSZ = [76 76];
for e = 1:numel(cidx1)
    rectangle('Position',[[cidx2(e) cidx1(e)]-BOXSZ/2 BOXSZ],'EdgeColor','g')
end
%% combine the hard threshold on the third score and the center finder to highlight the boxes we want
selectedBOX = imfill(~innerBOX_MASK,[cidx1 cidx2]) & innerBOX_MASK;
close all
tightenBoxFactor = 8;
imshow(selectedBOX,[])
selectedBOX = imerode(selectedBOX,strel('square',tightenBoxFactor));
R = regionprops(selectedBOX,'BoundingBox','Centroid');
imshow(bindVec(C2(:,:,1)),[]);
hold on
plot(cidx2,cidx1,'r*');
BOXSZ = [76 76];
for e = 1:numel(cidx1)
    rectangle('Position',[[cidx2(e) cidx1(e)]-BOXSZ/2 BOXSZ],'EdgeColor','g')
end
for e = 1:numel(R)
    rectangle('Position',R(e).BoundingBox,'EdgeColor','b')
end
%%
imshow(MASK,[]);