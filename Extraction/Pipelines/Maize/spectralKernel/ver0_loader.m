tPath = '/home/nate/Downloads/';
tFile = 'Oh43_shelled-1.hdr';
[wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo(tPath, tFile);
iFile = [tPath 'Oh43_shelled-1.raw'];
samplePart = readPart(iFile, spatial, spectral, framesUsed, 12);
%% read new data flour
tPath = '/mnt/tetra/nate/Jan182019/';
tFile = 'tril_1_20_1_Fluo.hdr';
[wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo(tPath, tFile);
iFile = [tPath 'tril_1_20_1_Fluo.raw'];
fid = fopen(iFile);
spectralImage_FL = readPart(fid, spatial, spectral, frames, 12);
fclose(fid);
SZ = size(spectralImage_FL);
tmp = reshape(spectralImage_FL,[prod(SZ(1:2)) SZ(3)]);
toComp1 = 5;
[S1,C1,U1,E1,L1,ERR1,LAM1] = PCA_FIT_FULL(double(tmp),toComp1);
C1 = reshape(C1,[SZ(1:2) toComp1]);
%% read new NIR
tPath = '/mnt/tetra/nate/Jan182019/';
tFile = 'tril_1_20_2_FX17.hdr';
[wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo(tPath, tFile);
iFile = [tPath 'tril_1_20_2_FX17.raw'];
fid = fopen(iFile);
spectralImage_NIR = readPart(fid, spatial, spectral, frames, 12);
fclose(fid);
spectralImage_NIR = flip(spectralImage_NIR,2);
SZ = size(spectralImage_NIR);
tmp = reshape(spectralImage_NIR,[prod(SZ(1:2)) SZ(3)]);
toComp2 = 5;
[S2,C2,U2,E2,L2,ERR2,LAM2] = PCA_FIT_FULL(double(tmp),toComp2);
C2 = reshape(C2,[SZ(1:2) toComp2]);
%%
[holeLocations_FL,R_seed_FL] = getHoleLocation(spectralImage_FL,true);
[holeLocations_NIR,R_seed_NIR] = getHoleLocation(spectralImage_NIR,true);
tform = fitgeotrans(holeLocations_NIR,holeLocations_FL,'projective');
%%
spectralImage_NIR_W = [];
SZ = size(spectralImage_FL);
R = imref2d(SZ(1:2));
for e = 1:size(spectralImage_NIR,3)
    spectralImage_NIR_W(:,:,e) = imwarp(spectralImage_NIR(:,:,e),tform,'OutputView',R);
    imshow(spectralImage_NIR_W(:,:,e),[]);
    drawnow
    e
end
C2_W = imwarp(C2,tform,'OutputView',R);
%%
close all
for e = 1:size(spectralImage_NIR,3)
    imshow(spectralImage_NIR(:,:,e)==255,[])
    drawnow
end
%%
reg1 = [];
reg2 = [];
close all
for e = 2:(numel(R_seed_FL)-1)

    ctmp1 = imcrop(C1(:,:,1:3),R_seed_FL(e).BoundingBox);
    ctmp2 = imcrop(C2_W(:,:,1:3),R_seed_FL(e).BoundingBox);

    D1 = [];
    for k = 1:size(spectralImage_NIR_W,3)
        D1(:,:,k) = imcrop(spectralImage_NIR_W(:,:,k),R_seed_FL(e).BoundingBox);
    end

    D2 = [];
    for k = 1:size(spectralImage_FL,3)
        D2(:,:,k) = imcrop(spectralImage_FL(:,:,k),R_seed_FL(e).BoundingBox);
    end

    uD1 = mean(D1,3);
    uD2 = mean(D2,3);

  
    MASK = bindVec(uD2) > .1;


    sz1 = size(D1);
    tmp1 = reshape(D1,[prod(sz1(1:2)) sz1(3)]);
    sz2 = size(D2);
    tmp2 = reshape(D2,[prod(sz2(1:2)) sz2(3)]);
    fidx = find(MASK == 1);


    reg1 = [reg1;tmp1(fidx,:)];
    reg2 = [reg2;tmp2(fidx,:)];


    Z = zeros(size(ctmp1));
    imshow([ctmp1;ctmp2;cat(3,bindVec(ctmp1(:,:,1)),Z(:,:,1),bindVec(ctmp2(:,:,1)))],[]);
    drawnow

    imshow(cat(3,bindVec(uD1),Z(:,:,1),bindVec(uD2)),[]);
    title(num2str(e))
    drawnow

end
%%
close all
[Xloadings,Yloadings,Xscores,Yscores, beta,pctVar,mse,stats,Weights] = plsregress(reg1,reg2,20);
preY = [ones(size(reg1,1),1) reg1]*beta;
for r = 1:10:size(preY,1)
    plot(preY(r,:),'r')
    hold on
    plot(reg2(r,:),'k');
    hold off
    drawnow
end
%%
uS_NIR_W = mean(spectralImage_NIR_W,3);
uS_FL = mean(spectralImage_FL,3);
Z = zeros(size(uS_FL));
imshow(cat(3,bindVec(uS_NIR_W),Z,bindVec(uS_FL)),[]);
%%
J1 = csvread('/mnt/tetra/macroPhor_corn/Histogram.csv');
J2 = csvread('/mnt/tetra/macroPhor_corn/NDVI.csv');
%% XLS read
label1 = {};
[D1 J1 R1] = xlsread('/mnt/tetra/macroPhor_corn/one.xlsx');
H1 = R1(1,:);
R1(1,:) = [];
for e = 1:size(R1,1)
    if ~isempty(R1(e,4))
        if R1{e,5} < 10
            spacer = '0'
        else
            spacer = '';
        end
        label1{e} = (['Corn ' num2str(R1{e,4}) spacer num2str(R1{e,5})]);
    end
    subD1 = cell2mat(R1(:,6:2:end));
end
label2 = {};
[D2 J2 R2] = xlsread('/mnt/tetra/macroPhor_corn/two.xlsx');
H2 = R2(1,:);
R2(1,:) = [];
for e = 1:size(R2,1)
    if ~isempty(R2(e,4))
        label2{e} = (['Corn ' num2str(R2{e,4}) num2str(R2{e,5})]);
    end
    subD2 = cell2mat(R2(:,6:2:end));
end
%%
FilePath = '/mnt/tetra/GWU/macroPhor_corn/';
FileList = {};
FileExt = {'hdr'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
close all
clear uS;
clear M;
parfor e = 1:numel(FileList)
    tic
    [p nm ext] = fileparts(FileList{e});
    rawFile = [p filesep nm '.raw'];
    [wavelengths, spatial, frames, spectral, tint, settings] = parseHdrInfo([p filesep],[nm ext]);
    rawFile = fopen(rawFile);
    samplePart = double(readPart(rawFile, spatial, spectral, frames, 12));
   
    [M{e} SP{e}] = makeKernelMask_ver0(samplePart);
    sz = size(samplePart);
    mini = reshape(samplePart,[prod(sz(1:2)) sz(3)]);
    uS(e,:) = mean(mini(find(M{e}),:),1);
    
    toc
    %imshow(M{e},[])
    %drawnow
    
    imgName{e} = nm;
    
end
%% fun regression
close all
X = [];
Y = [];
for e = 1:4
    Y = [Y;bindVec(SP{e}.loc(:,2))];
    %Y = [Y;(SP{e}.loc(:,1))];
    X = [X;SP{e}.data];
end
[XL,YL,XS,YS,BETA,PCTVAR] = plsregress(X,Y,3);
Yp = [ones(size(X,1),1) X]*BETA;
plot(Y(:),Yp(:),'.');
%% match
sidx = [];
for e = 1:numel(label1)
    sidx(e) = find(strcmp(label1{e},imgName));
end
imgName(sidx)
data = uS(sidx,:);
%data = data(randperm(50),:);
%%
for e = 1:numel(M)
    imshow(M{e},[])
    drawnow
    waitforbuttonpress
end
%% 
close all
figure;
[d1,p1] = corr(subD1,data);
mesh(d1);
view([0 90])
axis([1 size(data,2) 1 size(subD1,2)])
xlabel('wavelength')
ylabel('days')
title('NDVI')
colorbar


figure;
mesh(p1);
view([0 90])
axis([1 size(data,2) 1 size(subD1,2)])
xlabel('wavelength')
ylabel('days')
title('NDVI')
colorbar

figure;
[d2,p2] = corr(subD2,data);
mesh(d2);
view([0 90])
axis([1 size(data,2) 1 size(subD2,2)])
colorbar
xlabel('wavelength')
ylabel('days')
title('HISTOGRAM')

figure;
mesh(p2);
view([0 90])
axis([1 size(data,2) 1 size(subD1,2)])
xlabel('wavelength')
ylabel('days')
title('HISTOGRAM')
colorbar

%%

close all
figure;
[d1,p1] = corr(.5*(subD1(1:2:end,:) + subD1(2:2:end,:)),.5*(data(1:2:end,:) + data(2:2:end,:)));
mesh(d1);
view([0 90])
axis([1 size(data,2) 1 size(subD1,2)])
xlabel('wavelength')
ylabel('days')
title('NDVI')
colorbar


figure;
mesh(p1);
view([0 90])
axis([1 size(data,2) 1 size(subD1,2)])
xlabel('wavelength')
ylabel('days')
title('NDVI')
colorbar

figure;
[d1,p1] = corr(.5*(subD2(1:2:end,:) + subD2(2:2:end,:)),.5*(data(1:2:end,:) + data(2:2:end,:)));
mesh(d2);
view([0 90])
axis([1 size(data,2) 1 size(subD2,2)])
colorbar
xlabel('wavelength')
ylabel('days')
title('HISTOGRAM')

figure;
mesh(p2);
view([0 90])
axis([1 size(data,2) 1 size(subD1,2)])
xlabel('wavelength')
ylabel('days')
title('HISTOGRAM')
colorbar
%%
hypo = [];
close all
for hd = 1:25
    
    X1 = .5*(subD1(1:2:end,:) + subD1(2:2:end,:));
    X2 = .5*(subD2(1:2:end,:) + subD2(2:2:end,:));
    X = [X1 X2];
    X = X1;
    
    
    
  
    
    
    Y = .5*(data(1:2:end,:) + data(2:2:end,:));
    %Y = bsxfun(@times,Y,max(Y,[],2).^-1);  
    
    tmp = X;
    X = Y;
    Y = tmp;
    ind = setdiff(1:25,hd);
    trainX = X(ind,:);
    testX = X(hd,:);    
    trainY = Y(ind,:);
    testY = Y(hd,:);
    
    
    
    
    [XL,YL,XS,YS,BETA] = plsregress(trainX,trainY,3);
    
    
    hypo(hd,:) = [ones(1) testX]*BETA;
    
    
    
end

Y = .5*(data(1:2:end,:) + data(2:2:end,:));
Y = .5*(subD1(1:2:end,:) + subD1(2:2:end,:));
%Y = bsxfun(@times,Y,max(Y,[],2).^-1);  
[FUN,PFUN] = corr(hypo,Y);
csvwrite('/mnt/spaldingdata/Edgar/forGrantNOR.csv',[diag(FUN) diag(PFUN)]);
figure;
plot(diag(FUN))
figure;
plot(diag(PFUN))
figure;
mesh(FUN);
view([0 90])
axis([1 size(FUN,2) 1 size(FUN,1)])
colorbar

figure;
mesh(PFUN);
view([0 90])
axis([1 size(FUN,2) 1 size(FUN,1)])
colorbar

figure;
for e = 1:size(hypo,1)
    plot(hypo(e,:),'r')
    hold on
    plot(Y(e,:),'k')
    hold off
    waitforbuttonpress
end

figure;
plot(hypo','r')
hold on
plot(Y','k')




%% 
close all
figure;
mesh(corr(J1,uS));
view([0 90])
axis([1 size(uS,2) 1 size(J1,2)])
xlabel('wavelength')
ylabel('days')
title('Histogram')
colorbar
figure;
mesh(corr(J2,uS));
view([0 90])
axis([1 size(uS,2) 1 size(J1,2)])
colorbar
xlabel('wavelength')
ylabel('days')
title('NDVI')
%%
close all
figure
plot(uS')
legend(LEG)
%%
for e = 1:size(uS,1)
    uS(e,:) = bindVec(uS(e,:));
end
kidx = kmeans(uS(:,120:140),2)
plot(uS')

legend(LEG)
%%
for e = 1:size(samplePart,3)
    imshow(samplePart(:,:,e),[]);
    drawnow
end
%%
close all
uS = mean(samplePart,2);
imshow(squeeze(uS),[])
%%
uS = mean(samplePart,2);
[S C U E L ERR LAM] = PCA_FIT_FULL(squeeze(uS),3);
%%
sz = size(samplePart);
F = reshape(samplePart,[prod(sz(1:2)) sz(3)]);
[T] = PCA_REPROJ(F,E,U);
T = reshape(T,[sz(1:2) 3]);
%%
for k = 1:size(T,3)
    T(:,:,k) = bindVec(T(:,:,k));
end
imshow(T,[]);
%%
close all
sz = size(samplePart);
F = reshape(samplePart,[prod(sz(1:2)) sz(3)]);
[T] = PCA_REPROJ(F,E,U);
BK = kmeans(T,3);
T = reshape(T,[sz(1:2) 3]);
BK = reshape(BK,[sz(1:2)]);
rgbBK = label2rgb(BK);
imshow(rgbBK,[]);
%%



    
