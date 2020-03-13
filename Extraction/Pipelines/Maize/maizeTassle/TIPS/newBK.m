
%% try maxwell
I = double(imread('~/Downloads/DSC_0081.JPG'))/255;
Io = double(I)/255;
reSZ = 10;
BOX = findTasselCropBox(Io,reSZ,.5);
boxImg = imcrop(Io, BOX);
boxMask = thresholdTasselImage(boxImg, reSZ);
finalMask = tasselMask2OriginalSize(boxMask, BOX, Io);



I = imcrop(I,[]);
%% find white rectangle
reSZ = 10;
%%
close all
reSZ = 10;
for e = 3:numel(FileList)

    if ~contains(FileList{e},'background')

        Io = double(I)/255;
        %I = permute(I,[2 1 3]);
        %I = flip(I,1);


        BOX = findTasselCropBox(Io,reSZ,.5);
        I = imcrop(Io,BOX);

        Lab = rgb2lab(I);
        %Lab = rgb2gray(I);


        isColor = mean(abs(Lab(:,:,2:3)),3);
        %isColor = max(abs(Lab(:,:,2:3)),[],3);
        raw = imcomplement(isColor);
        %raw2 = imcomplement(Lab(:,:,1));
        %raw3 = mean(cat(3,raw,raw2),3);
        %raw = raw3;
        
        
        
        toProcess = raw;
        oldSZ = size(toProcess);
        toProcess = imresize(toProcess,1/reSZ);
        toProcess = imfilter(toProcess,fspecial('gaussian',[31 31],7),'replicate');
        toProcess = imdilate(toProcess,strel('line',101,0));
        toProcess = imfilter(toProcess,fspecial('gaussian',[31 31],7),'replicate');
        toProcess = imresize(toProcess,oldSZ);
       
        
      

        tassel = toProcess - raw;
        %tassel = imcomplement(raw);
        tassel = bindVec(tassel);
        tasselM = tassel > graythresh(tassel);
        tasselM = imclose(tasselM,strel('disk',5,0));
        
        tasselM = imclearborder(tasselM);
        
        tasselM = bwlarge(tasselM);
        tasselM = imclose(tasselM,strel('disk',5,0));
        
        tidx = [];
        [tidx(:,2),tidx(:,1)] = find(tasselM);
       
        if ~isempty(tidx)
            tidx(:,1) = tidx(:,1) + round(BOX(1));
            tidx(:,2) = tidx(:,2) + round(BOX(2));
        end
        T = sub2ind([size(Io,1) size(Io,2)],tidx(:,2),tidx(:,1));
        tasselM = zeros(size(Io,1),size(Io,2));
        tasselM(T) = 1;
        out = flattenMaskOverlay(Io,logical(tasselM));
        %tasselM = imfill(tasselM,'holes');
        %tasselM = imopen(tasselM,strel('disk',7,0));

   

        imshow(out,[]);
        hold on
        rectangle('Position',BOX,'EdgeColor','r')
        drawnow
    end
end
