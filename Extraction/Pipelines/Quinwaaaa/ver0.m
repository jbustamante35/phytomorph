I = imread('/home/nate/IMG_6632.JPG');

FilePath = '/mnt/tetra/nate/QU/test/AllPhotos/';
FileList = {};
FileExt = {'tif','TIF','JPG'};
FileList = gdig(FilePath,FileList,FileExt,1);
oPath = '/mnt/tetra/nate/QU/test/return/';
mkdir(oPath);
K = [];
for e = 1:numel(FileList)
    I = imread(FileList{e});
    if e == 1
        K = I;
    else
        I = imresize(I,[size(K,1) size(K,2)]);
        K = (.9*K+.1*I);
    end
    imshow(K,[]);
    drawnow
end
%{
    Lab = rgb2lab(I);
    HSV = rgb2hsv(I);
    isColor = mean(abs(Lab(:,:,2:3)),3);
    isColorS = std((Lab(:,:,2:3)),1,3);
    isColor = isColor.*Lab(:,:,1);
    %isColor = bindVec(isColor);
    MASK = isColor > 1000;%graythresh(isColor);
    %MASK = imclose(MASK,strel('disk',31,0));
    %MASK = imclearborder(MASK);
    MASK = bwareaopen(MASK,50000);
    %MASK = bwlarge(MASK,3);
    fidx = find(MASK);
    R = regionprops(MASK,'MajorAxisLength','MinorAxisLength','PixelIdxList','BoundingBox','Orientation','Centroid','Perimeter','EulerNumber','Area');
    for r = 1:numel(R)
        for k = 1:size(I,3)
            tmp = I(:,:,k);
            uRGB(r,k) = mean(tmp(R(r).PixelIdxList));
        end
    end

    %
    close all
    out = flattenMaskOverlay(I,MASK,.4,'r');
    imshow(out,[]);
    hold on
    TH = linspace(-pi,pi,200);

    for r = 1:numel(R)
        CIR = [];
        R(r).Orientation = R(r).Orientation*pi/180;
        CIR(:,1) = R(r).MajorAxisLength/2*cos(TH);
        CIR(:,2) = R(r).MinorAxisLength/2*sin(TH);
        ROT = [[cos(R(r).Orientation) sin(R(r).Orientation)];[-sin(R(r).Orientation) cos(R(r).Orientation)]];
        CIR = (ROT*CIR')';
        rectangle('Position',R(r).BoundingBox,'EdgeColor','r');
        rectangle('Position',[R(r).BoundingBox(1:2) 100 100],'FaceColor',uRGB(r,:)/255);
        plot(CIR(:,1)+R(r).Centroid(1),CIR(:,2)+R(r).Centroid(2),'c');
    end
    drawnow
    [pth,nm,ext] = fileparts(FileList{e});
    saveas(gca,[oPath nm '.jpg']);
    
    
    
    
    
end
%}
%%
close all
out = flattenMaskOverlay(I,MASK);
imshow(out,[])
%%