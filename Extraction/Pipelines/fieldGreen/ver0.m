FilePath = '/mnt/tetra/nate/forCraine/RGB/';
rgbFileList = {};
FileExt = {'JPG'};
rgbFileList = gdig(FilePath,rgbFileList,FileExt,1);
%%
%%
FilePath = '/mnt/tetra/nate/forCraine/NIR1/';
nirFileList = {};
FileExt = {'JPG'};
nirFileList = sdig(FilePath,nirFileList,FileExt,1);
%% sample for the rgb
CSNIR = [];
for e = 1:numel(rgbFileList)
    I = imread(nirFileList{e});
    Lab = rgb2lab(I);
    
    I = double(I);
    
    if exist('GMM_NIR')
        sz = size(I);
        tmp = double(reshape(Lab,[prod(sz(1:2)) sz(3)]));
        kidx = GMM_NIR.cluster(tmp);
        kidx = reshape(kidx,sz(1:2));
        plantMask = kidx == 2;
        out = flattenMaskOverlay(double(I)/255,plantMask);
       
        
        
        DISKMask = Lab(:,:,1) > 98;
        DISKMask = bwlarge(DISKMask);
        DISKMask = imerode(DISKMask,strel('disk',2,0));
        
        for k = 1:3
            tmp = I(:,:,k);
            urgb(k) = mean(tmp(find(DISKMask)));
        end
        delta1 = 255 - urgb(1);
        I(:,:,1) = I(:,:,1) + delta1;
        
        
        NDVI = (I(:,:,1) - I(:,:,2)).*(I(:,:,1) + I(:,:,2)).^-1;
        
        uNDVI = mean(NDVI(find(plantMask)));
        sNDVI = std(NDVI(find(plantMask)));
        
        out = flattenMaskOverlay(out,DISKMask,.5,'b');
        
        imshow(out,[]);
        waitforbuttonpress
    end
    
    
    
    
    
    tmp = imresize(Lab,.25);
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);
    tmp = tmp(randperm(size(tmp,1)),:);
    CSNIR = [CSNIR;tmp(1:5000,:)];
    imshow(I,[]);
    drawnow
    e
end
%%
options = statset('Display','iter');
GMM_NIR = fitgmdist(double(CSNIR),3,'Options',options);
%%
close all
for e = 1:1
    singleGreenFieldImage(rgbFileList{e},GMM,'./output/')
end
%%
func = @(X,T)singleGreenFieldImage(X,GMM,GMM_NIR,T,'./output/');
funcP = partialFunction(func,'fieldImage_WSU');
funcP.publish();
%%
h1 = figure;
h2 = figure;
h3 = figure;
FL = [];
PER = [];
for e = 1:numel(rgbFileList{1})
    tmpPER = [];
    uNDVI = [];
    for t = 1:numel(rgbFileList)
        I = imread(rgbFileList{t}{e});
        FL = [FL I];
        
        [tmpPER(t),uNDVI(t)] = func(nirFileList{t}{e},'NIR');
        figure(h1)
        plot(tmpPER,'r')
        figure(h2)
        plot(uNDVI,'g')
        drawnow
        e
        figure(h3);
        imshow(I,[]);
        
    end
    PER(t,:) = tmpPER;
end
%%

func(nirFileList{1},'NIR');
%% sample for the rgb
CS = [];
for e = 1:numel(rgbFileList)
    I = imread(rgbFileList{e});
    Lab = rgb2lab(I);
    
    
    
    imshow(I,[]);
    if exist('GMM')
        sz = size(I);
        tmp = double(reshape(Lab,[prod(sz(1:2)) sz(3)]));
        kidx = GMM.cluster(tmp);
        kidx = reshape(kidx,sz(1:2));
        plantMask = kidx == 2;
        out = flattenMaskOverlay(double(I)/255,plantMask);
       
        
        
        DISKMask = Lab(:,:,1) > 99;
        DISKMask = bwlarge(DISKMask);
        out = flattenMaskOverlay(out,DISKMask,.5,'b');
        
        imshow(out,[]);
        waitforbuttonpress
    end
    
    
    
    
    
    tmp = imresize(Lab,.25);
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);
    tmp = tmp(randperm(size(tmp,1)),:);
    CS = [CS;tmp(1:5000,:)];
    imshow(I,[]);
    drawnow
    e
end
%%
options = statset('Display','iter');
GMM = fitgmdist(double(CS),3,'Options',options);




