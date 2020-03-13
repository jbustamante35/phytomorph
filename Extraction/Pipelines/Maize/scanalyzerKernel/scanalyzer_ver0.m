%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local process of remote data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather PNG files from lucia user
dataPath = ['/iplant/home/lucia/corn_seeds%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.png')
        KerFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mass download data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target = '/mnt/spaldingdata/nate/octerineDataStore/argeliaDataSets/';
source = '/iplant/home/lucia/corn_seeds_color_2019%';
imassDownload(source,target,'png',12)
%% create FinFileList
for e = 1:numel(KerFileList)
    [pth{e},nm,ext] = fileparts(KerFileList{e});
end
UQ = unique(pth);
for e = 1:numel(UQ)
    kp(e) = contains(UQ{e},'vis');
end
UQ = UQ(find(kp));
for e = 1:numel(UQ)
    FinFileList{e} = {};
end
cnt =1;
for e = 1:numel(KerFileList)
    for u = 1:numel(UQ)
        if contains(KerFileList{e},UQ{u})
            FinFileList{u}{end+1} = KerFileList{e};
            ScanFileList{cnt} = KerFileList{e};
            cnt = cnt + 1;
        end
    end
end
%% try all local
for e = 1:numel(FinFileList)
    scanalyzerMain(FinFileList{e},tform,tform2,'/mnt/tetra/nate/Argelia_model/',[]);
end

%%
%%%%%%%%%%%%%%%%%%%%%%
%%
FilePath = '/mnt/tetra/nate/cornScan/corn_seeds_images_08-07-2018/';
FileList = {};
FileExt = {'tif','TIF','png'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%













IDX = 1:2:numel(FileList2);
IDX2 = 2:2:numel(FileList2);
for n = 1:numel(FileList)
    iF = FileList{n};
    fF = FileList2{IDX(n)};
    nF = FileList2{IDX2(n)};

    scanalyzerMain(iF,fF,nF,tform,tform2,'./output/',[]);
end
%%
FilePath = '/mnt/tetra/nate/Argelia_model/';
FilePath = '/mnt/tetra/nate/lucia/analyses/';
tFileList = {};
FileExt = {'tif'};
tFileList = gdig(FilePath,tFileList,FileExt,1);
FilePath = '/mnt/tetra/nate/Argelia_model/';
FilePath = '/mnt/tetra/nate/lucia/analyses/';
jFileList = {};
FileExt = {'json'};
jFileList = gdig(FilePath,jFileList,FileExt,1);
FilePath = '/mnt/tetra/nate/Argelia_model/';
FilePath = '/mnt/tetra/nate/lucia/analyses/';
mFileList = {};
FileExt = {'mat'};
mFileList = gdig(FilePath,mFileList,FileExt,1);
 %% intersect against mat files
matFileListCheck = {};
cnt = 1;
RGBIDX = [];
for e = 1:numel(tFileList)
    if contains(tFileList{e},'RGB')
        matFileListCheck{cnt} = strrep(strrep(strrep(tFileList{e},[filesep 'RGB'],''),'SQUARE','extra'),'.tif','.mat');
        RGBIDX(cnt) = e;
        cnt = cnt + 1;
    end
end
[~,IDX,~] = intersect(matFileListCheck,mFileList);
tFileList = tFileList(RGBIDX(IDX));
%%
cnt = 1;
SZ = [];
S = [];
S2 = [];
S3 = [];
CS = [];
K = [];
rnd = 800;
close all
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
h5 = figure;
EMB_M = [];
KER_M = [];

for e = 1:numel(tFileList)
    e
    if contains(tFileList{e},'RGB')
        try
            I = double(imread(tFileList{e}))/255;
            I2 = double(imread(strrep(tFileList{e},'RGB','FLUO')))/255;
            I3 = double(imread(strrep(tFileList{e},'RGB','IR')))/255;
            jFN = strrep(strrep(tFileList{e},[filesep 'RGB'],''),'SQUARE','jdoc');
            mFN = strrep(strrep(tFileList{e},[filesep 'RGB'],''),'SQUARE','extra');
            jFN = strrep(jFN,'.tif','.json');
            mFN = strrep(mFN,'.tif','.mat');
            iFN = strrep(tFileList{e},'/mnt/tetra/nate/Argelia_model/','/iplant/home/lucia/corn_seeds/');

            didx = strfind(tFileList{e},'__');
            iFN_N = tFileList{e}(didx(1)+2:end);
            iFN_N(end-10:end) = [];
            iFN_N = [iFN_N '.png'];



            rm = strfind(iFN,'RGB');
            iFN(rm:end) = [];

            %{
            key = strrep(iFN,'/iplant/home/lucia/corn_seeds/','');
            kidx1 = find(contains(ScanFileList,key));
            kidx2 = find(contains(ScanFileList,iFN_N));
            KK = intersect(kidx1,kidx2);
            iFN = ScanFileList{KK};
            J = imread(iFN);
            %}
            J = ones(size(J));
            objD = load(mFN);
            tmpD = loadjson(jFN);
            fI = imfilter(I,fspecial('gaussian',[21 21],7),'replicate');
            I = imresize(I,[300 790]);
            fI = imresize(fI,[300 790]);
            I2 = imresize(I2,[300 790]);
            I3 = imresize(I3,[300 790]);
            S(:,:,:,cnt) = I;
            S2(:,:,:,cnt) = I2;
            S3(:,:,:,cnt) = I3;
            CON(:,:,cnt) = tmpD.kernelDoc.interpolatedContour.data;


            drawnow
            SZ(cnt,:) = size(I);



            tmpcs = reshape(fI,[300*790 3]);
            kidx = kernelGMM.cluster(tmpcs);
            kidx = reshape(kidx,[300 790]);
            %kidx = ones([300 790]);
            K(:,:,cnt) = kidx;


            EMB = K(:,:,cnt)==1;
            EMB = bwlarge(EMB);


            eB = bwboundaries(EMB);
            eB = eB{1};
            eBi = [];
            eBi(:,1) = ba_interp2(objD.G(:,:,1),eB(:,2),eB(:,1));
            eBi(:,2) = ba_interp2(objD.G(:,:,2),eB(:,2),eB(:,1));

            kBi = tmpD.kernelDoc.interpolatedContour.data;




            figure(h1);
            imshow(J,[]);
            hold on
            plot(kBi(:,2),kBi(:,1),'r');
            hold on
            plot(eBi(:,2),eBi(:,1),'b');
            CM = mean(kBi,1);
            plot(CM(2),CM(1),'k*');
            plot(kBi(1,2),kBi(1,1),'k*');
            TAN = kBi(1,:) - CM;
            TAN = TAN / norm(TAN);
            NOR = [-TAN(2) TAN(1)];
            quiver(CM(2),CM(1),TAN(2),TAN(1),1000);
            quiver(CM(2),CM(1),-TAN(2),-TAN(1),1000);
            quiver(CM(2),CM(1),NOR(2),NOR(1),1000);
            quiver(CM(2),CM(1),-NOR(2),-NOR(1),1000);
            hold off

            kMASK = poly2mask(kBi(:,2),kBi(:,1),size(J,1),size(J,2));

            % measure embryo shape
            iKER = [];
            iKERi = [];
            [iKER(:,1),iKER(:,2)] = find(kMASK);
            iKERi(:,1) = ba_interp2(objD.G(:,:,1),iKER(:,2),iKER(:,1));
            iKERi(:,2) = ba_interp2(objD.G(:,:,2),iKER(:,2),iKER(:,1));
            iKERi = bsxfun(@minus,iKER,(kBi(1,:)));
            Fx = round(iKERi*[(-TAN)' (NOR)']);
            Fx = unique(Fx,'rows');
            [~,sidx] = sort(Fx(:,1));
            Fx = Fx(sidx,:);
            UQ = unique(Fx(:,1));
            KER_shape = [];
            for i = 1:numel(UQ)
                KER_shape(i,:) = [UQ(i) sum(Fx(:,1)==UQ(i))];
            end

            % measure embryo shape
            iEMB = [];
            iEMBi = [];
            [iEMB(:,1),iEMB(:,2)] = find(EMB);
            iEMBi(:,1) = ba_interp2(objD.G(:,:,1),iEMB(:,2),iEMB(:,1));
            iEMBi(:,2) = ba_interp2(objD.G(:,:,2),iEMB(:,2),iEMB(:,1));
            iEMBi = bsxfun(@minus,iEMBi,(kBi(1,:)));
            Fx = round(iEMBi*[(-TAN)' (NOR)']);
            Fx = unique(Fx,'rows');
            [~,sidx] = sort(Fx(:,1));
            Fx = Fx(sidx,:);
            UQ = unique(Fx(:,1));
            EMB_shape = [];
            for i = 1:numel(UQ)
                EMB_shape(i,:) = [UQ(i) sum(Fx(:,1)==UQ(i))];
            end

            DELTA = [];
            for i = 1:size(EMB_shape,1)
                xPOS = EMB_shape(i,1);

                fidx = find(KER_shape(:,1) == xPOS);
                if ~isempty(fidx)
                    DELTA(i,:) = [xPOS EMB_shape(i,2) - KER_shape(fidx,2)];
                end
            end

            fidx = find(abs(DELTA(:,2)) <= max(7,min(abs(DELTA(:,2)))));
            EMB_shape(1:fidx(end),:) = [];
            fidx = find(KER_shape(:,1) == EMB_shape(1,1));
            EMB_shape = [KER_shape(1:fidx - 1,:);EMB_shape];

            YV = linspace(0,EMB_shape(end,1),1000);
            XV = interp1(EMB_shape(:,1),EMB_shape(:,2),YV,'linear','extrap');

            YW = linspace(0,KER_shape(end,1),1000);
            XW = interp1(KER_shape(:,1),KER_shape(:,2),YW,'linear','extrap');

            KER_M = cat(3,KER_M,[YW' XW']);
            EMB_M = cat(3,EMB_M,[YV' XV']);


            if cnt > 10
                kerLEN = squeeze(KER_M(end,1,:));
                embLEN = squeeze(EMB_M(end,1,:));
                X = squeeze(KER_M(:,2,:))';
                Y = squeeze(EMB_M(:,2,:))';



                % stack matrix for analysis
                Z = [X';Y';kerLEN';embLEN']';
                % find outlier and remove
                [zS,zC,zU,zE,zL,zERR,zLAM] = PCA_FIT_FULL(Z,3);
                TF = isoutlier(zC);
                Z(any(TF,2),:) = [];
                X(any(TF,2),:) = [];
                Y(any(TF,2),:) = [];
                kerLEN(any(TF,2),:) = [];
                embLEN(any(TF,2),:) = [];
                % move forward with clean dataset

                [zS,zC,zU,zE,zL,zERR,zLAM] = PCA_FIT_FULL(Z,100);
                [sweepD] = sweepPCA(zC,zE,zU,diag(zLAM).^.5,1:size(zE,2),5);

                comKer = zE(1:1000,:);
                comEmb = zE(1001:2000,:);

                [comQKer,R] = qr(comKer,0);
                [comQEmb,R] = qr(comEmb,0);


                indX = X - (comQKer*(comQKer'*X'))';
                indY = Y - (comQEmb*(comQEmb'*Y'))';


                [wxS,wxC,wxU,wxE,wxL,wxERR,wxLAM] = PCA_FIT_FULL(X,5);
                [wyS,wyC,wyU,wyE,wyL,wyERR,wyLAM] = PCA_FIT_FULL(Y,5);

                [xS,xC,xU,xE,xL,xERR,xLAM] = PCA_FIT_FULL(indX,5);
                [yS,yC,yU,yE,yL,yERR,yLAM] = PCA_FIT_FULL(indY,5);

                swX = sweepPCA(xC,xE,wxU,diag(xLAM).^.5,1,10);
                swY = sweepPCA(yC,yE,wyU,diag(yLAM).^.5,1,10);



                figure(h5)
                for para = 1:2
                    figure(h5)
                    subplot(1,2,para)
                    sweepZ = squeeze(sweepD(para,:,:))';
                    nX = sweepZ(1:1000,:);
                    nY = sweepZ(1001:2000,:);
                    lX = sweepZ(2001,:);
                    lY = sweepZ(2002,:);
                    for p = 1:size(nX,2)
                        lXi = linspace(0,lX(p),1000);
                        lYi = linspace(0,lY(p),1000);

                        plot(lXi,nX(:,p),'k')
                        hold on
                        plot(lXi,-nX(:,p),'k')
                        plot(lYi,nY(:,p),'r')
                        plot(lYi,-nY(:,p),'r')
                    end
                    drawnow
                    hold off
                end
                %{
                for para = 1:3
                    figure
                    sweepZ = squeeze(sweepD(para,:,:))';
                    nX = sweepZ(1:1000,:);
                    nY = sweepZ(1001:2000,:);
                    lX = sweepZ(2001,:);
                    lY = sweepZ(2002,:);




                    for p = 1:size(nX,2)
                        lXi = linspace(0,lX(p),1000);
                        lYi = linspace(0,lY(p),1000);

                        plot(lXi,nX(:,p),'k')
                        hold on
                        plot(lXi,-nX(:,p),'k')
                        plot(lYi,nY(:,p),'r')
                        plot(lYi,-nY(:,p),'r')
                    end
                    title(['PCA' num2str(para)])
                    drawnow
                    hold off
                axis([0 850 -600 600])
                end
                %}

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [xS,xC,xU,xE,xL,xERR,xLAM] = PCA_FIT_FULL(X,3);
                [yS,yC,yU,yE,yL,yERR,yLAM] = PCA_FIT_FULL(Y,3);
                [A,B,r,U,V,stats] = canoncorr(xC,yC);

                toProbe = 1;
                [fS,fC,fU,fE,fL,fERR,fLAM] = PCA_FIT_FULL([U(:,toProbe),V(:,toProbe)],1);
                newC = linspace(-fLAM,fLAM,5);
                tmpData = PCA_BKPROJ(newC',fE,fU);
                tmpDataU = [tmpData(:,1) zeros(size(tmpData,1),2)];
                tmpDataV = [tmpData(:,2) zeros(size(tmpData,1),2)];
                simX = tmpDataU*inv(A);
                simY = tmpDataV*inv(B);
                simX = PCA_BKPROJ(simX,xE,xU);
                simY = PCA_BKPROJ(simY,yE,yU);
                embX = linspace(0,mean(embLEN),1000);
                kerX = linspace(0,mean(kerLEN),1000);
                figure(h4);
                subplot(1,2,1)
                plot(embX',simY','r');
                hold on
                plot(embX',-simY','r');
                plot(kerX',simX','k');
                plot(kerX',-simX','k');
                hold off
                title([num2str(corr(fC,squeeze(kerLEN))) '--' num2str(corr(fC,squeeze(embLEN)))])
                kerLEN = squeeze(KER_M(end,1,:));
                embLEN = squeeze(EMB_M(end,1,:));
                kerLEN(any(TF,2),:) = [];
                embLEN(any(TF,2),:) = [];


                toProbe = 2;
                [fS,fC,fU,fE,fL,fERR,fLAM] = PCA_FIT_FULL([U(:,toProbe),V(:,toProbe)],1);
                newC = linspace(-fLAM,fLAM,5);
                tmpData = PCA_BKPROJ(newC',fE,fU);
                tmpDataU = [zeros(size(tmpData,1),1) tmpData(:,1) zeros(size(tmpData,1),1)];
                tmpDataV = [zeros(size(tmpData,1),1) tmpData(:,2) zeros(size(tmpData,1),1)];
                simX = tmpDataU*inv(A);
                simY = tmpDataV*inv(B);
                simX = PCA_BKPROJ(simX,xE,xU);
                simY = PCA_BKPROJ(simY,yE,yU);
                embX = linspace(0,mean(embLEN),1000);
                kerX = linspace(0,mean(kerLEN),1000);
                figure(h4);
                subplot(1,2,2)
                plot(embX',simY','r');
                hold on
                plot(embX',-simY','r');
                plot(kerX',simX','k');
                plot(kerX',-simX','k');
                hold off
                title([num2str(corr(fC,squeeze(kerLEN))) '--' num2str(corr(fC,squeeze(embLEN)))])

                toProbe = 2;
                [fS,fC,fU,fE,fL,fERR,fLAM] = PCA_FIT_FULL([U(:,toProbe),V(:,toProbe)],1);
                newC = linspace(-fLAM,fLAM,5);
                tmpData = PCA_BKPROJ(newC',fE,fU);
                tmpDataU = [zeros(size(tmpData,1),1) tmpData(:,1) zeros(size(tmpData,1),1)];
                tmpDataV = [zeros(size(tmpData,1),1) tmpData(:,2) zeros(size(tmpData,1),1)];
                simX = tmpDataU*inv(A);
                simY = tmpDataV*inv(B);
                simX = PCA_BKPROJ(simX,xE,xU);
                simY = PCA_BKPROJ(simY,yE,yU);
                embX = linspace(0,mean(embLEN),1000);
                kerX = linspace(0,mean(kerLEN),1000);
                figure(h4);
                subplot(1,2,2)
                plot(embX',simY','r');
                hold on
                plot(embX',-simY','r');
                plot(kerX',simX','k');
                plot(kerX',-simX','k');
                hold off
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end





            figure(h3)
            plot(KER_shape(:,1),KER_shape(:,2))
            hold on
            plot(EMB_shape(:,1),EMB_shape(:,2))
            hold off

            EmbryoArea = size(Fx,1);

            IDX = randperm(size(tmpcs,1));
            CS = [CS;tmpcs(IDX(1:rnd),:)];
            SNG = cat(1,mean(I,4),mean(I2,4),mean(I3,4),repmat(K(:,:,cnt)==2,[1 1 3]));
            MUL = cat(1,mean(S,4),mean(S2,4),mean(S3,4),repmat(mode(K,3)/NG,[1 1 3]));

            figure(h2);
            imshow([SNG MUL],[]);

            figure(h5);
            drawnow
            pause(1)
            cnt = cnt + 1;
        catch ME
            getReport(ME)
        end
    end
end
%%
commonBasis
%% sample the color space first
CS = [];
for e = 1:numel(tFileList)
    if contains(tFileList{e},'RGB')
        I = double(imread(tFileList{e}))/255;
        I = imresize(I,[300 790]);
        IDX = randperm(size(tmpcs,1));
        tmpcs = reshape(I,[300*790 3]);
        CS = [CS;tmpcs(IDX(1:rnd),:)];
        e
    end
end
%% cluster before running the embryo code above
NG = 3;
options = statset('Display','iter');
kernelGMM = fitgmdist(CS,NG,'Options',options);
%%
S = permute(S,[1 2 4 3]);
S2 = permute(S2,[1 2 4 3]);
S3 = permute(S3,[1 2 4 3]);
sz1 = size(S);
s1 = reshape(S,[prod(sz1(1:3)) sz1(4)]);
s2 = reshape(S2,[prod(sz1(1:3)) sz1(4)]);
s3 = reshape(S3,[prod(sz1(1:3)) sz1(4)]);
%%
close all
SKIP = 100;
plot(s1(1:SKIP:end,3),s3(1:SKIP:end,1),'.')
%%
sS = [];
for e = 1:size(S,4)
    sS(:,:,:,e) = imresize(S(:,:,:,e),.1);
end
sS(:,:,:,any(TF,2)) = [];
SSS = size(sS);
sS = reshape(sS,[prod(SSS(1:3)) SSS(4)]);
[Se Ce Ue Ee Le ERRe LAMe] = PCA_FIT_FULL_T(sS,3);
%%
close all
SD = diag(LAMe).^.5;
[sweepD] = sweepPCA(Ce',Ee,Ue',SD(1:2),[1 2],15);
for pc = 1:size(sweepD,1)
    for s = 1:size(sweepD,2)
        data = sweepD(pc,s,:);
        data = reshape(data,SSS(1:3));
        data = imresize(data,4);
        imshow(data,[]);
        title(num2str(pc))
        drawnow
    end
end
%%
close all
h1 = figure;
h2 = figure;
DD = [];
[zS,zC,zU,zE,zL,zERR,zLAM] = PCA_FIT_FULL(Z,15);
for h = 1:10
    for k = 1:15
        [Train, Test] = crossvalind('HoldOut', size(U,1), .2);
        %trainX = [X(Train,:) Y(Train,:)];
        %testX = [X(Test,:) Y(Test,:)];
        trainX = [zC(Train,:)];
        testX = [zC(Test,:)];
        trainY = sS(:,Train)';
        testY = sS(:,Test)';
        TestIDX = find(Test);
        

        [Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,Weights] = plsregress(trainX,trainY,k);


        YI = [ones(size(testX,1),1) testX]*beta;

        DELTA = YI - sS(:,Test)';
        DD(h,k) = mean(sum(DELTA.^2,2));
        YI = reshape(YI',[SSS(1:3) sum(Test)]);

        figure(h1);
        errorbar(mean(DD,1),std(DD,1,1))
        

        figure(h2);

        for q = 1:10%size(YI,3)
            imshow(cat(1,imresize(YI(:,:,:,q),4),imresize(reshape(sS(:,TestIDX(q)),SSS(1:3)),4)),[]);
            title(num2str(k))
            %waitforbuttonpress
            pause(.1)
        end
        

        figure(h1);
        drawnow
    end
end
%%
FilePath = '/mnt/tetra/nate/cornScan/corn_seeds_images_FLUO_Aand_NIR_08_12_18/';
FileList2 = {};
FileExt = {'tif','TIF','png'};
FileList2 = sdig(FilePath,FileList2,FileExt,1);
%%
iF = FileList{1};
%% for first tform
I1 = imread(FileList{1}{1});
I2 = imread(FileList2{1}{1});
cpselect(I1,I2)
tform = fitgeotrans(fixedPoints2,movingPoints2,'projective');
R = imref2d(size(I2));
B = imwarp(I2,R,tform);
%%
I1 = imread(FileList{1}{1});
I2 = imread(FileList2{2}{1});
cpselect(I1,I2)
tform2 = fitgeotrans(fixedPoints3,movingPoints3,'projective');
%% look at sweeps of kernel alone
D1.d = X;
D2.d = Y;
D3.d = sS';
D4.d = [kerLEN embLEN];
B={D1 D2 D3 D4};
options=PCAGCA(B);
options.ExpVarLim=5; %each component has to explain at least 5% of the variation in the blocks
options.Rlim=0.50; %the canonical correlation has to be above 0.9 in order to define a component as common
%options.autoselect=1;
model=PCAGCA(B,options);

