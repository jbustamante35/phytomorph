% orginal reads from guosheng images
I = imread('/home/nate/Downloads/p1_mica_20160823_index_ndvi.tif');
%I = imread('/mnt/spaldingdata/Drone_Imagery/Arlington_2018/2018_07_23/2018_07_23.tif');
I = imread('/home/nate/GridNoGrid.tif');
%% read(s) for Dayane's images
I = imread('/home/nate/orthomosaic.tif');
imshow(I,[]);
%%
I = I(:,:,1:3);
%% 
imshow(I,[]);
%%
G = rgb2gray(I);
%%
imshow(G,[]);
%% rotate manual for now
close all
rI = imrotate(I,1);
imshow(rI,[]);
G = rgb2gray(rI);
%% major crop level one
close all
[subRange,box] = imcrop(G);
subColor = imcrop(rI,box);
close all
imshow(subRange,[]);
%% refine crop box level two
close all
view1 = [400 600];
sub1 = subColor(1:view1(1),1:view1(2),:);
[cp(1,1),cp(1,2),v1] = impixel(sub1);
sub2 = subColor((end-view(1)+1):end,1:view1(2),:);
[cp(2,1),cp(2,2),v2] = impixel(sub2);
offset1 = size(subColor,1)-view1(1);
cp(2,2) = cp(2,2) + offset1;
sub3 = subColor(1:view1(1),(end-view(2)+1):end,:);
[cp(3,1),cp(3,2),v3] = impixel(sub3);
offset2 = size(subColor,2)-view1(2);
cp(3,1) = cp(3,1) + offset2;
sub4 = subColor((end-view(1)+1):end,(end-view(2)+1):end,:);
[cp(4,1),cp(4,2),v4] = impixel(sub4);
cp(4,1) = cp(4,1) + offset2;
cp(4,2) = cp(4,2) + offset1;
%% show clicks
close all
imshow(subColor,[]);
hold on
plot([cp(1,1) cp(2,1)],[cp(1,2) cp(2,2)],'r')
plot([cp(3,1) cp(4,1)],[cp(3,2) cp(4,2)],'r')
plot([cp(1,1) cp(3,1)],[cp(1,2) cp(3,2)],'b')
plot([cp(2,1) cp(4,1)],[cp(2,2) cp(4,2)],'b')

%% perform rotation and show rotated clicks on image
close all

angleFix(1) = atan2(cp(2,1)-cp(1,1),cp(2,2)-cp(1,2))*180/pi;
angleFix(2) = atan2(cp(4,1)-cp(3,1),cp(4,2)-cp(3,2))*180/pi;
angleFix(3) = atan2(cp(3,1)-cp(1,1),cp(3,2)-cp(1,2))*180/pi - 90;
angleFix(4) = atan2(cp(4,1)-cp(2,1),cp(4,2)-cp(2,2))*180/pi - 90;
uA = mean(angleFix);


%uA = 90;
testColor = imrotate(subColor,uA,'crop');
testGray = imrotate(subRange,uA,'crop');
ruA = uA*pi/180;
M = [[cos(ruA);sin(ruA)],[-sin(ruA);cos(ruA)]];
tcp = bsxfun(@minus,cp,flip(hsize(subRange),2));
tcp = (M'*tcp')';
tcp = bsxfun(@plus,tcp,flip(hsize(testGray),2));
imshow(testColor,[]);
hold on

plot([tcp(1,1) tcp(2,1)],[tcp(1,2) tcp(2,2)],'r')
plot([tcp(3,1) tcp(4,1)],[tcp(3,2) tcp(4,2)],'r')
plot([tcp(1,1) tcp(3,1)],[tcp(1,2) tcp(3,2)],'b')
plot([tcp(2,1) tcp(4,1)],[tcp(2,2) tcp(4,2)],'b')
pM = tcp;
tmp = pM(3:4,:);
pM(3:4,:) = flip(tmp,1);
fieldMask = poly2mask(pM(:,1),pM(:,2),size(testGray,1),size(testGray,2));
imshow(fieldMask,[]);
R = regionprops(fieldMask);
subRange = imcrop(subRange,R.BoundingBox);
subColor = imcrop(subColor,R.BoundingBox);
%%
close all
imshow(subColor,[])
%%
figure;
imshow(subColor,[])
%% find the range
% the threshold values were hard coded for guoshengs data
% dilate size should be the ~ the length of the row
dilateValue = 800; % guoshengs data
dilateValue = 200; % dayane
thresholdValue  = -80; % guosheng data
thresholdValue = -120;
filterSize = 100; % for guosheng
filterSize = 20; % for Dayane

close all
sig = mean(subRange,2);
figure;
plot(sig,'b')
sig = imfilter(sig,ones(filterSize,1)/filterSize,'symmetric');
hold on;
plot(sig,'r')
rangeD = imdilate(sig,ones(dilateValue,1)) == sig & sig > thresholdValue;
plot(rangeD*max(sig),'k');
waitforbuttonpress


raIDX = find(rangeD);
plot(-100*rangeD,'k')
close all
figure;
imshow(subRange,[]);
hold on
plot((-30*sig)-min(-30*sig),1:size(sig,1),'r')
plot((5000*rangeD),1:size(sig,1),'y')
%% display the range crop lines
%axis([0 size(rangeD,1) 0 2]);
%{
IM = imreconstruct(sig-100,sig);
sig = sig - IM;
rangeD = sig > 0;
%}
raIDX = find(rangeD);
rangeD = repmat(rangeD>0,[1 size(subRange,2)]);
rangeD = imdilate(rangeD,strel('disk',11,0));
out = flattenMaskOverlay(subRange,rangeD);
figure;
imshow(rangeD,[]);
figure;
imshow(out,[]);
%imwrite(out,'/mnt/tetra/nate/for_IOWA_CORN/large.tiff');
%%
subI = double(subI)/255;
%%
%%
close all
sig = fft(subI-mean(subI,2),[],2);
sig2 = mean(abs(sig),1);
ssig2 = imfilter(sig2,fspecial('average',[1 21]),'replicate');
figure;
plot(ssig2,'k')
hold on
plot(sig2);
figure;
sig3 = sig2 - ssig2;
plot(sig3);
sig4 = imdilate(sig3,strel('disk',51,0)) == sig3;
cutoff = 50;
sig4(1:cutoff) = 0;
fidx = find(sig4);
hold on
plot(sig4*30);
%%
BW = roipoly(I);
%%

close all
subI = imcrop(I);
subI = bindVec(subI);
imshow(double(subI),[]);
%%
%%
sig = [];
for r = 1:(numel(raIDX)-1)
    subI = double(subRange(raIDX(r):raIDX(r+1),:))/255;
    close all
    [f,phase] = findMajorFrequency(subI,50,140^-1,[],0,[],[]);
    
end
%%
F = mean(sig,1);
FF = std(sig,1,1);
%%
close all
figure;
imshow(subRange,[],'Border','tight');
hold on
rangeD = imdilate(sig,ones(800,1)) == sig & sig > -80;
%plot((-30*sig)-min(-30*sig),1:size(sig,1),'r')
plot((5000*rangeD),1:size(sig,1),'y')

for r = 1:(numel(raIDX)-1)
    for e = 1:numel(RW{r})
        plot(RW{r}(e).Centroid(1)*ones(size(subI{r},1),1),raIDX(r)+(1:size(subI{r},1)),'r');
        hold on
    end
    %{
    for e = 1:numel(MR{r})
        plot(MR{r}(e).Centroid(1)*ones(size(subI{r},1),1),raIDX(r)+(1:size(subI{r},1)),'y');
        hold on
    end
    %}
end
%% row finding
%for r rows
subI = {};
RW = {};
for r = 1:(numel(raIDX)-1)
    % crop out the "plot"
    subI{r} = imresize(double(subRange(raIDX(r):raIDX(r+1),:))/255,1);
    subC{r} = imresize(double(subColor(raIDX(r):raIDX(r+1),:,:))/255,1);
    close all
    
    [f,phase] = findMajorFrequency_ver2(subI{r},[4 30]);
    
    %[f,phase] = findMajorFrequency(subI{r},50,(.5*140*2)^-1,0); guoshend
    %data
    
    %[f,phase] = findMajorFrequency(subI,50,50^-1,0);
    [mask] = generateRowMask(subI{r},f,phase);
    out = flattenMaskOverlay(double(bindVec(subI{r})),mask>.8);
    mask2 = repmat(mean(mask,1),[size(mask,1) 1]);
    
    out2 = flattenMaskOverlay(double(bindVec(subI{r})),bindVec(mask2)<.02);
    out2 = flattenMaskOverlay(subC{r},bindVec(mask2)<.02);
    
    imshow(out,[])
    title([num2str(f^-1) '--' num2str(round(size(subI{r},2)*f))]);
    waitforbuttonpress
    
    subO = out2(:,1:600,:);
    close all
    figure;
    imshow(subO,[]);
    
    
    
    waitforbuttonpress
    
    %{
    RW{r} = regionprops(mask2 < -.2,'centroid');
    hold on
    for e = 1:numel(RW{r})
        plot(RW{r}(e).Centroid(1)*ones(size(subI{r},1),1),1:size(subI{r},1),'r');
        hold on
    end
    
    %}
    
    %imwrite(out2,['/mnt/tetra/nate/for_IOWA_CORN/small' num2str(r) '.tiff']);
    %title(num2str(f^-1))
    %figure;
    %plot(mean(subI{r},1));
    %waitforbuttonpress
end
%%

%%
%v = VideoWriter('/mnt/tetra/nate/FFT2.avi');
%v.Quality = 100;
%open(v)
%close all
blockNUM = 25;
FINAL_BLOCK_SIZE = 30;
for row = 1:numel(subI)
    BIG = zeros(size(subI{row}));
    Tphase = [];
    for blockNUM = 5:5:FINAL_BLOCK_SIZE
        WIN = round((1/f)*blockNUM);
        for s = 1:10:(size(subI{row},2)-WIN)

            tempy = subI{row}(:,s:(s+WIN));
            IDXf = round(size(tempy,2)*f)+1;
            tic
            [Tf,Tphase(s,:)] = findMajorFrequency(tempy,50,140^-1,0);
            toc
            [Tmask] = generateRowMask(tempy,f,Tphase(s,:));
            
            Tmask = repmat(mean(Tmask,1),[size(Tmask,1) 1]);



            tic
            BLOCK = BIG(:,s:(s+WIN));
            BLOCK = mean(cat(3,BLOCK,Tmask),3);
            BIG(:,s:(s+WIN)) = BLOCK;
            toc



            iBOOL = imerode(mean(BIG,1),strel('disk',31)) == mean(BIG);
            
            R = regionprops(iBOOL,'centroid');

            
            if blockNUM < 2* FINAL_BLOCK_SIZE
                if mod(s,100) == 1
                    imshow(subI{row},[],'Border','tight');
                    hold on
                    for r = 1:numel(R)
                        %plot(R(r).Centroid(1)*ones(size(subI{row},1),1),linspace(1,(size(subI{row},1)),size(subI{row},1)),'y');
                        hold on
                    end


                    %plot(s*ones(size(subI{row},1),1),1:size(subI{row},1),'g');
                    %plot((s+WIN)*ones(size(subI{row},1),1),1:size(subI{row},1),'r');
                    %{
                    out = flattenMaskOverlay(subI,BIG > .6);
                    imshow(out,[]);
                    %}

                    plot(-mean(BIG,1)*100 + size(subI{row},1)-100,'m','LineWidth',2)



                     for e = 1:numel(RW{row})
                        %plot(RW{row}(e).Centroid(1)*ones(size(subI{row},1),1),linspace(1,size(subI{row},1),size(subI{row},1)),'c');
                        hold on
                    end
                    drawnow
                    hold off

                    
                    %frame = getframe(gcf);
                    %writeVideo(v,frame);
                    
                end
            end
            s
            (size(subI{row},2)-WIN)
            MR{row} = R;
            %{
            imshow(BIG,[]);
            drawnow

            %}
            %{
            imshow(tempy,[]);
            drawnow
            %}
        end
    end
    %waitforbuttonpress
end
%close(v)
%%
fBIG = imfilter(BIG,fspecial('gaussian',[31 31],12));

%%
close all
    out = flattenMaskOverlay(subI,BIG > .1);
    imshow(out,[]);
%%
close all
for e = 1:1
    

    out2 = flattenMaskOverlay(double(bindVec(subI)),bindVec(mask2)<.02);
    imshow(out2,[]);
    hold on
    sig = -mean(BIG,1);
    IM = imreconstruct(sig-.1,sig);
    sig = sig - IM;
    sig = repmat(sig > 0,[size(BIG,1) 1]);
    waitforbuttonpress
    out2 = flattenMaskOverlay(out2,sig,.5,'b');
    imshow(out2,[]);
    drawnow
    %plot(mean(BIG,1)*500+500,'r')
    waitforbuttonpress
    
    %plot(mean(mask2,1)*500+500,'b')
end
%%
ridx = find(any(BW,2));
PAD = 51;
phun = zeros(size(BW));
for r = PAD:(numel(ridx)-PAD)
    
    
    tic
    strip = I((ridx(r)-PAD):(ridx(r)+PAD),:);
    Mstrip = BW((ridx(r)-PAD):(ridx(r)-PAD),:);
    midx = find(all(Mstrip,1));
    subI = strip(:,midx);
    toc
    
    
    
    [~,BOX] = imcrop(I);
  
    
    tic
    [block] = getAblock(I,BOX);
    toc
    
    block = imfilter(block,fspecial('gaussian',[15 15],2),'replicate');
    
    tic
    [f,phase] = findMajorFrequency(block,50,40^-1,0);
    toc
    
    tic
    [mask] = generateRowMask(block,f,phase);
    toc
    
    rowMask = mask < -.85;
    imshow(rowMask,[]);
    out = flattenMaskOverlay(double(bindVec(block)),rowMask,.5,'r');
    figure;
    imshow(out,[]);
    
    
    sig1 = mean(block,1);
    sig1 = zscore(sig1);
    
    sig2 = mean(rowMask,1);
    sig2 = zscore(sig2);
    figure
    plot(sig1,'g')
    hold on
    plot(sig2,'r')
    
    figure;
    imshow(out(:,1:300,:),[]);
    
    
    
    
    
    phun(ridx(r),midx) = mask((end-1)/2,:);
    
    %imshow(phun,[]);
    %drawnow
    r
    numel(ridx)
end
%%
close all
wow = mean(I.*BW,2);
plot(wow.*(wow>.4))
widx = find(wow > .2);
wf = fft(wow(widx));
plot(abs(wf))
bf = 45/size(wf,1);
bm = cos(2*pi*(1:size(I,1))*bf);
bm = repmat(bm',[1 size(I,2)]);
imshow(bm > .95,[]);
out = flattenMaskOverlay(double(I.*BW),bm>.95,.5,'b');
bwf = bm > .95;
%%
phunM = phun < -.95;
%%
out = flattenMaskOverlay(double(bindVec(BW.*I)),logical(BW.*(bwf|phunM)),.7,'r');
%%
close all
plot(mean(subI,1))
%%
close all
phase = angle(sig(:,113));
for r = 1:size(subI,1)
    what(r,:) = cos(2*pi*(1:size(subI,2))*T^-1 + phase(r));
end
imshow(what,[])
msk = what > .8;
out = flattenMaskOverlay(double(subI),msk);
imshow(out,[]);

