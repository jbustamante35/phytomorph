%%

D = readtext('/home/nate/Downloads/SD_2016_genotypeNames.csv');
%D = readtext('/home/nate/Downloads/sd_2017_names.csv');
%%
FilePath = '/mnt/tetra/nate/stomataTopoData/Accessions_2016/';
FileList = {};
FileExt = {'nms'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc
%%
[S C U E L ERR LAM] = PCA_FIT_FULL_T(A,3);
%%
for k = 1:3
    iC(:,:,k) = bindVec(col2im(C(k,:),[51 51],[512 512]));
end
close all
imshow(iC,[]);
%%
X = table;
cnt = 1;
Y = [];
nSZ = [512 512];
Xi = [];
SS = {};
IDX = {};
for e = 1:161%numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    ridx = find(strcmp(D(:,2),strrep(['16' nm],'','')));
    %X(cnt,'filename') = {FileList{e}};
    %Xi(:,:,cnt) = imresize(imread(FileList{e}),nSZ);
    %A = im2colF(imresize(imread(FileList{e}),nSZ),[51 51],[1 1]);
    
    %SS{e} = mean(A,2);
    %IDX{e} = ridx;
    
    if ~isempty(ridx)
        X(cnt,'filename') = {FileList{e}};
        %Xi(:,:,cnt) = imresize(imread(FileList{e}),nSZ);
        A = im2colF(Xi(:,:,cnt),[51 51],[1 1]);
        SS(:,cnt) = sum(A,2);
        X(cnt,'Count') = {D{ridx,3}};
        Y(cnt) = D{ridx,3};
        cnt
        cnt = cnt + 1;
        
    end
    
    e
end
%%
XS = [];
for e = 1:size(X,1)
    XS(:,:,1,e) = imread(X.filename{e});
    e
end
%%
XS = [];
YS = [];
for e = 1:numel(IDX)
    if ~isempty(IDX{e})
        YS = [YS D{IDX{e}(1),3}];
        XS = [XS SS{e}];
    end
end
%%
Xi = reshape(Xi,[size(Xi,1) size(Xi,2) 1 size(Xi,3)]);
[Yi,U,S] = zscore(YS);
%%
%cnt = 1;
close all
figure
for e = 161:numel(FileList)
    
    [pth,nm,ext] = fileparts(FileList{e});
    ridx = find(strcmp(D(:,2),strrep(['16' nm],'','')));
    
    
    
    if ~isempty(ridx)
        
        
        
        I = imread(FileList{e});
        tmp = trainedNet2.predict(I);
        tmp2 = trainedNetNN.predict(I);
        tmp3 = D{ridx,3};
        
       
        imshow(I,[]);
        title([num2str(tmp) ',' num2str(tmp2) '-->' num2str(tmp3)]);
        waitforbuttonpress
        tmp = (inputdlg('k'));
        handY(cnt) = str2num(tmp{1});
        handX(:,:,1,cnt) = I;
        
        
        
        
        cnt = cnt + 1;
    end
    
    e
     
end
%%
[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(XS',YS',31);
Ypre = beta'*[ones(1,size(XS,2));XS];
%%
I = imread(FileList{1});
SZ = 1;

layers = [imageInputLayer([512 512],'Normalization','none');
          convolution2dLayer([23,45],1);
          reluLayer
          maxPooling2dLayer(21,'Stride',20)
          fullyConnectedLayer(1)
          regressionLayer];
%{ 
layers = [imageInputLayer([512 512],'Normalization','none');
      convolution2dLayer([23,45],1);
      fullyConnectedLayer(1)
      regressionLayer];
          %}
  
%{
layers = [imageInputLayer([512 512],'Normalization','none');
          convolution2dLayer([23,45],1);
          reluLayer
          maxPooling2dLayer(11,'Stride',5);
          convolution2dLayer([11,11],1);
          reluLayer
          maxPooling2dLayer(11,'Stride',5)
          fullyConnectedLayer(1)
          regressionLayer];
%}
      
%{
layers = [imageInputLayer([512 512],'Normalization','none');
          convolution2dLayer([23,45],1);
          reluLayer
          maxPooling2dLayer(51,'Stride',23);
          convolution2dLayer([11,11],9);
          reluLayer
          maxPooling2dLayer(3,'Stride',5)
          fullyConnectedLayer(1)
          regressionLayer];
%}
%{
[q1 q2 V] = impixel(XS(:,:,2));
W(:,:,1,1) = XS(q2-11:q2+11,q1-22:q1+22,2);
[q1 q2 V] = impixel(XS(:,:,6));
W(:,:,1,2) = XS(q2-11:q2+11,q1-22:q1+22,6);
[q1 q2 V] = impixel(XS(:,:,100));
W(:,:,1,3) = XS(q2-11:q2+11,q1-22:q1+22,100);
%}
for r = 1:10
    [q1 q2 V] = impixel(XS(:,:,r+1));
    W(:,:,1,r) = XS(q2-11:q2+11,q1-22:q1+22,r+1);
end
%%
for l = 1:size(W,4)
    tmp = W(:,:,1,l);
    
    [d1 d2] = gradient(tmp);
    tmp = (d1.^2 + d2.^2).^.5;
    
    
    tmp = tmp - mean(tmp(:));
    s = std(tmp(:));
    tmp = tmp / s;
    Wg(:,:,1,l) = tmp;
end
%%
for h = 1:size(handX,4)
    tmp = handX(:,:,1,h);
    
    [d1 d2] = gradient(tmp);
    tmp = (d1.^2 + d2.^2).^.5;
    
    
    handXg(:,:,1,h) = tmp;
end
%%
for h = 1:size(handX,4)
    tmp = handX(:,:,1,h);
    Ui(h) = mean(tmp(:));
    tmp = tmp - mean(tmp(:));
    Si(h) = std(tmp(:));
    tmp = tmp / std(tmp(:));
    
    
    handXn(:,:,1,h) = tmp;
end
%%

layers = [imageInputLayer([512 512 1],'Normalization','none');
          convolution2dLayer([23,45],1);
          reluLayer
          maxPooling2dLayer(11,'Stride',5)
          
          convolution2dLayer([11,11],9);
          reluLayer;
          maxPooling2dLayer(11,'Stride',5);
          
          fullyConnectedLayer(1)
          regressionLayer];
      
      
[handYi,U,S] = zscore(handY);
close all
%layers(2).Weights = cat(3,mean(W,4),mean(Wg,4),zeros(size(W,1),size(W,2)));
%layers(2).Weights = mean(W,4);
%layers(2).Weights = cat(3,trainedNetRaw.Layers(2).Weights,trainedNetGrad.Layers(2).Weights,zeros(size(trainedNetGrad.Layers(2).Weights)));
layers(2).Weights = mean(W,4);
%layers(2).Weights = W;
layers(2).Bias = rand(1,1,1)*.0001 + 1;
%layers(2).Bias = cat(3,trainedNetRaw.Layers(2).Bias,trainedNetGrad.Layers(2).Bias,zeros(size(trainedNetGrad.Layers(2).Bias)));
%trainX = cat(3,handX(:,:,:,1:end-20),handXg(:,:,:,1:end-20),zeros(size(handXg(:,:,:,1:end-20))));
trainX = handXn(:,:,:,1:end-20);
trainY = handY(1:end-20);
%testX = cat(3,handX(:,:,:,end-19:end),handXg(:,:,:,end-19:end),zeros(size(handXg(:,:,:,end-19:end))));
testX = handXn(:,:,:,end-19:end);
testY = handY(end-19:end);
options = trainingOptions('sgdm','Shuffle','every-epoch','ValidationFrequency',10,'ValidationData',{testX,testY'},'MiniBatchSize',25,'InitialLearnRate',0.00000001,'Plots','training-progress','MaxEpochs',200,'ExecutionEnvironment','parallel');
%AUG = imageDataAugmenter('RandXScale',[1/SZ 1/SZ],'RandYScale',[1/SZ 1/SZ]);
%imagesource = augmentedImageSource(size(I)/SZ,X,'DataAugmentation',AUG);
%trainedNet = trainNetwork(imagesource,layers,options);
%trainedNetNN = trainNetwork(XS,Y',layers,options);

trainedNetNNN = trainNetwork(trainX,trainY',layers,options);
%%
for e = 1:4
    trainedNetNNN = trainNetwork(trainX,trainY',trainedNetNNN.Layers,options);
    e
end

%%
close all
PPY = S*trainedNetNNi.predict(handX)+U;
PPY2 = trainedNetNNN.predict(handX);
plot(handY,PPY,'.')
hold on
plot(0:200,'r-')
figure
plot(handY,PPY2,'.');
hold on
plot(0:200,'r-')
figure;
plot(handY,mean([PPY,PPY2],2),'.')
hold on
plot(0:200,'r-')
%%
options = trainingOptions('sgdm','MiniBatchSize',128,'InitialLearnRate',0.000000001,'Plots','training-progress','MaxEpochs',45,'ExecutionEnvironment','parallel');
trainedNet2 = trainNetwork(XS,Y',trainedNet2.Layers,options);
%%
options = trainingOptions('sgdm','MiniBatchSize',64,'InitialLearnRate',0.001,'Plots','training-progress','MaxEpochs',45,'ExecutionEnvironment','parallel');
trainedNet3 = trainNetwork(newX,pY',trainedNet2.Layers,options);
%% try network
close all
Y = [];
pY = [];
pY2 = [];
pY3 = [];
pYM = [];
newX = [];
cnt = 1;
parfor e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    ridx = find(strcmp(D(:,2),strrep(['16' nm],'','')));
    
    tic
    if ~isempty(ridx)
        I = imread(FileList{e});
        [d1,d2] = gradient(I);
        dI = (d1.^2 + d2.^2).^.5;
        nI = (I - mean(I(:))) / std(I(:));
        %I = imresize(I,nSZ);
        %A = im2colF(I,[51 51],[1 1]);
        %A = sum(A,2);
        %A = reshape(A,[51 51]);
        tmp = trainedNetGrad.predict(dI);
        tmp2 = trainedNetRaw.predict(I);
        tmp3 = trainedNetNor.predict(nI);
        
        %X(cnt,'filename') = {FileList{e}};
        %X(cnt,'Count') = {D{ridx,3}};
        
        
        Y(e) = D{ridx,3};
        
        pY(e) = tmp;
        pY2(e) = tmp2;
        pY3(e) = tmp3;
        
        
        %pYM = mean([pY',pY2',pY3'],2);
        IDX{e} = ridx;
        %{
        %if abs(pY(cnt) - Y(cnt)) > 10
        %if cnt == 38
            th = figure;
            imshow(I,[]);
            hold on
            
            title([num2str(pY(cnt)) ',' num2str(pY2(cnt)) ',' num2str(pY3(cnt)) ':' num2str(pYM(cnt)) '-->' num2str(Y(cnt))])
            waitforbuttonpress
            pause(1);
            close(th)
        %end
        %}
        %}
        %newX(:,:,:,cnt) = I;
        %{
        hold on
        plot(Y,pY,'b.')
        plot(Y,pY2,'c.')
        plot(Y,pY3,'k.')
        plot(Y,pYM,'r.')
        hold on
        plot(0:200,'r--')
        hold off
        drawnow
        cnt = cnt + 1;
        %}
    else
        Y(e) = NaN;
        
        pY(e) = NaN;
        pY2(e) = NaN;
        pY3(e) = NaN;
        IDX{e} = ridx;
    end
    toc*numel(FileList)/12/60
    
    
    
end
%%
DK = D;
%%
pYM = nanmean([pY',pY2',pY3'],2);
for e = 1:numel(IDX)
    idx = IDX{e};
    if ~isempty(idx)
        for r = 1:numel(idx)
            D{idx(r),7} = pY(e);
            D{idx(r),8} = pY2(e);
            D{idx(r),9} = pY3(e);
            D{idx(r),10} = pYM(e);
        end
    end
end
D{1,7} = 'new1';
D{1,8} = 'new2';
D{1,9} = 'new3';
D{1,10} = 'new4';
%%
close all
plot(Y,pY3,'.')
%%
close all
plot(Y,pYM,'.')
%%  remove Q1
for e = 1:size(D,1)
    if contains(D{e,11},'Q1')
        rm(e) = true;
    else
        rm(e) = false;
    end
end
D(rm,:) = [];
%%
%genoData = D(2:end,end);
genoData = D(2:end,10);
genoNames = unique(genoData);
%phenoData = cell2mat(D(2:end,4));
phenoData = cell2mat(D(2:end,7));
sampleData = D(2:end,4);
%%
ksdensity(phenoData)
%%
close all
gm = {};
AIC = [];
OUTPUT = {};
for e = 1:numel(genoNames)
    fidx = find(strcmp(genoData,genoNames{e}));
    subP = phenoData(fidx,:);
    options = statset('Display','off');
    samp = sampleData(fidx);
    
    for g = 1:2
        flag = 0;
        cnt = 1;
        while flag == 0
            try
                %isoutlier(subP,'method','gesd')
                
                gm{e,g} = fitgmdist(subP,g,'Options',options,'Replicates',10,'RegularizationValue',.01);
                AIC(e,g) = gm{e,g}.AIC;
                %AIC(e,g) = gm{e,g}.BIC;
                %AIC(e,g) = -gm{e,g}.NegativeLogLikelihood;
                flag = 1;
            catch
                flag = 0;
                cnt = cnt + 1;
                if cnt > 20
                    flag = 1;
                    AIC(e,g) = inf;
                end
                
            end
        end
        
    end
    
    
    [Y,X] = ksdensity(subP);
    
    %plot(X,Y,'b');
    %hold on
    [~,midx] = min(AIC(e,:));
    A = gm{e,midx};
    kidx = A.cluster(subP);
   
    [~,MIDX] = max(A.ComponentProportion);
    k = kidx == MIDX;
    
    CL = {'k'};
    for s = 1:numel(A.ComponentProportion)
        Yi = normpdf(X,A.mu(s),A.Sigma(s).^.5);
     %   plot(X,Yi,CL{1});
        if s == MIDX
       %    plot(X,Yi,'r'); 
            valueP = A.mu(s);
        end
    end
    
    %title([genoNames{e} '-->' num2str(e)])
    %pause(.01);
    %waitforbuttonpress
    
    %close all
    %{
    ksdensity()
    drawnow
    hold on
    %}
    
    OUTPUT{e,1} = genoNames{e};
    OUTPUT{e,2} = valueP;
    OUTPUT{e,3} = mean(subP);
    e
    numel(genoNames)
end
%%
close all
good = find(AIC(:,1) < AIC(:,2));
good = min(AIC,[],2);
for e = 1:numel(good)
    fidx = find(strcmp(genoData,genoNames{good(e)}));
    subP = phenoData(fidx,:);
   
    
    ksdensity(subP)
    drawnow
    hold on
    waitforbuttonpress
end
