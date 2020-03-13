%D = readtext('~/k_1.ubh','\t');
D = csvread('~/collection-070918_131242_vert_lidar');
%P = cell2mat(D(:,2:end));
P = D(300:1000,2:end);
P(P >= (65533 - 10)) = 0;
%%
TH = linspace(-3*pi/4,3*pi/4,size(P,2));
H = bsxfun(@times,P,cos(TH));
W = bsxfun(@times,P,sin(TH));
TM = repmat((1:size(P,1))',[1 size(H,2)]);
%%
figure;
imshow(P,[]);
%%
kidx = W > 0 & W < 800;
W = W(kidx);
H = H(kidx);
TM = TM(kidx);
%%
close all
plot3(W(:),H(:),TM(:),'.')
axis equal

%%
close all
ksdensity(W(:))
%%
close all
mag = 1;
[iSZ,offset] = getISZ(H,W,TM,mag);
%iSZ(1) = 2*iSZ(1);
WTH = [0 800];
%WTH = [0 300];
[R,offset] = slicer(H,W,TM,WTH,iSZ,offset,mag);
[R,GR] = groundRemove(R);

WTH = [800 1600];
%WTH = [200 400];
[G] = slicer(H,W,TM,WTH,iSZ,offset,mag);
[G] = groundRemove(G,GR);

WTH = [1600 inf];
%WTH = [600 inf];
[B] = slicer(H,W,TM,WTH,iSZ,offset,mag);
[B] = groundRemove(B,GR);

%B = zeros(size(R));
RGB = cat(3,R,G,B);
RGB = imresize(RGB,[size(RGB,1) size(RGB,1)]);
figure;

vRGB = imfilter(RGB,fspecial('disk',11),'replicate');
for k = 1:size(RGB,3)
    vRGB(:,:,k) = bindVec(vRGB(:,:,k));
end
imshow(flip(vRGB,1),[]);
figure;
imshow(flip(RGB(:,:,1) >0 ,1),[]);
%%
s = imcrop(flip(vRGB,1));
imshow(s,[]);
hold on
pH1 = sum(s(:,:,1),2);
plot(pH1, 1:numel(pH1),'r');
%%
close all
plot(sum(RGB(:,:,2),1))
%%
close all
HH = sum(RGB(:,:,2)>0,2);
plot(HH)

%%
slice = flip(RGB(100:200,:,1),1);
sliceF = im2col(slice,[size(slice,1) 21],'sliding');
[sS,sC,sU,sE,sL,sERR,sLAM] = PCA_FIT_FULL_T(sliceF,3);
sC = padarray(sC,[0 10],0,'both');


slice = flip(RGB(300:600,:,1),1);
sliceF = im2col(slice,[size(slice,1) 21],'sliding');
[sS,sC2,sU,sE,sL,sERR,sLAM] = PCA_FIT_FULL_T(sliceF,3);
sC2 = padarray(sC2,[0 10],0,'both');
%%
figure;
plot(sum(slice(:,:,1),1)) 
%%
close all
figure;
imshow(RGB,[]);
%%
figure;
imshow(flip(RGB,1),[]);
hold on
mag = 10;
pause(1)
plot(mag*sC(1,:) + size(RGB,1)*95/100,'c');
plot(mag*sC2(1,:) + size(RGB,1)*95/100,'m');
sig = sC2(1,:)<10;
sig = sum(RGB(200:end,:,1),1);
sig = imfilter(sig,fspecial('disk',11));
sig = sig > 8;

sig = ~bwareaopen(~sig,30);
plot(mag*sig + size(RGB,1)*95/100,'y');
R = regionprops(logical(sig),'PixelIdxList','BoundingBox');
IMG = [];
for r = 1:numel(R)
    
    BOX = R(r).BoundingBox;
    BOX(4) = size(RGB,1);
    rectangle('Position',BOX,'EdgeColor','r');
    tmp = imcrop(RGB(:,:,1),BOX);
    
    tmp = imresize(tmp(1:200,:),[200 500]);
    imshow(tmp);
    drawnow
    IMG = cat(3,IMG,tmp);
end



%%
close all
rpH = round(pH);
rpW = round(pW);
rpT = round(pTM);
rpH = rpH - min(rpH) + 1;
rpW = rpW - min(rpW) + 1;
rpT = rpT - min(rpT) + 1;
max(rpH)
max(rpW)
max(rpT)

I = zeros(max(rpH),max(rpT));
IDX = sub2ind(size(I),rpH,rpT);
I(IDX) = rpW;

figure;
imshow(I,[]);
%%
close all
I = imresize(I,[size(I,1) size(I,1)]);
imshow(flip(I,1),[])

figure;
imshow(imfilter(flip(I,1),fspecial('disk',[11]),'replicate'),[]);
%%

for j = 1:size(I,2)
    tmp = I(:,j) > 0;
    tmp = [ones(size(tmp)) tmp ones(size(tmp))];
    tmp = imfill(logical(tmp),[1,2],8) - tmp;
    tmp = tmp(:,2);
    sidx = find(tmp);
    if ~isempty(sidx)
        I(:,j) = circshift(I(:,j),-sidx(end));
    end
    
end

figure;
imshow(imfilter(flip(I,1),fspecial('disk',[11]),'replicate'),[]);
%%
imshow(I(end-200:end,:),[]);
%%
figure;
imshow(I(1:200,:),[]);
figure;
imshow(I,[])
figure;
imshow(I==0,[])
J = I == 0;
figure;
for j = 1:size(J,2)
    tmp = I(:,j) > 0;
    tmp = [ones(size(tmp)) tmp ones(size(tmp))];
    tmp = imfill(logical(tmp),[1,2],8) - tmp;
    tmp = tmp(:,2);
    sidx = find(tmp);
    if ~isempty(sidx)
        I(:,j) = circshift(I(:,j),-sidx(end));
    end
    
end
figure;imshow(I,[]);
figure;
imshow(imfilter(I,fspecial('gaussian',[31 31]),'replicate'),[]);
%%
I = imresize(I,[size(I,1) size(I,1)]);
%%
close all
imagesc(bindVec(flip(I,1)))
colormap('jet')
%%
figure;
imshow(I,[]);
%%
close all
sig = I(100:300,100:end);
F = abs(fft(bsxfun(@minus,sig,mean(sig,2)),[],2));
plot(mean(F,1));
%%


%%
close all
for t = 1:size(I,2)
    plot(I(end-200:end,t))
    drawnow
    
end



