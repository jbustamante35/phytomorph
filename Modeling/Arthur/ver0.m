% load  data
%M = bfopen('/mnt/tetra/nate/arthur/Experiment-983.czi');
%M = bfopen('/mnt/tetra/nate/arthur/Experiment-974.czi');
M = bfopen('/mnt/tetra/nate/arthur/Experiment-974.czi');
%% register leaves
J = bfopen('/mnt/tetra/nate/arthur/Experiment-974.czi',1,1,1,1);
%% read tiff stack
info = imfinfo('/mnt/tetra/nate/arthur/041719_1uMflg22_spl001 (1).tif');
for e = 1:numel(info)
    M(:,:,:,e) = imread('/mnt/tetra/nate/arthur/041719_1uMflg22_spl001 (1).tif',e);
    e
end
%% convert to gray scale STACK
close all
STACK = [];
for e = 1:size(M,4)
    
    tmp = single(rgb2gray(M(:,:,:,e)))/255;
    if e == 1
        [tmp,BOX] = imcrop(tmp);
    else
        tmp = imcrop(tmp,BOX);
    end
    
    tmp = imfilter(tmp,fspecial('gaussian',[41 41],11),'replicate');
    STACK(:,:,e) = tmp;
    %{
    imshow(tmp,[]);
    drawnow
    %}
    e
end
%% conver to single and 255
STACK = single(STACK)/255;
%% STACK data
close all
h1 = figure;
h2 = figure;
cnt = 1;
STACK = zeros([size(M{1}{1,1}) size(M{1},1)],'single');
for t = 1:1:size(M{1},1)
    tmp = single(M{1}{t,1});
    tmp(tmp == max(tmp(:))) = min(tmp(:));
    tmp = imfilter(tmp,fspecial('gaussian',[41 41],11),'replicate');
    %{
    %tmp = bindVec(tmp);
    min(tmp(:))
    max(tmp(:))
    figure(h1);
    %subplot(1,2,1);
    imshow(tmp,[]);
    title(num2str(t))
    %}
    STACK(:,:,cnt) = tmp;
    cnt = cnt + 1;
    %{
    figure(h1);
    subplot(1,2,2);
    imshow(std(STACK,1,3),[]);
    
    
    drawnow
    %}
    t
    size(M{1},1)
    
end
%%
%{
%% remove zero frames - for down sample
STACK(:,:,cnt:end) = [];
%% play movie again
close all
for e = 100:20:size(STACK,3)
    imshow(STACK(:,:,e),[]);
    title(num2str(e))
    drawnow
end
%% for sub sample
STACK = STACK(:,:,100:1000);
%% diff
dS = diff(STACK,1,3);
%%
close all
hist(dS(:),linspace(min(dS(:)),max(dS(:)),1000));
%% grad
THRESH = 3;
R = [min(dS(:)) max(dS(:))];
close all
for e = 1:size(dS,3)
    sig = dS(:,:,e);
    sig = imfilter(sig,fspecial('gaussian',[21 21],11),'replicate');
    %Tsig = adaptthresh(sig,1,'ForegroundPolarity','dark');
    %MSK = imbinarize(sig,Tsig);
    MSK = sig  > THRESH;
    MSK = MSK .* pM;
    out = flattenMaskOverlay(pM.*bindVec(STACK(:,:,e)),logical(MSK));
    %imshow(sig,R);
    imshow(out,[]);
    %contour(sig)
    drawnow
end
%}
clear M
%% make mask for plant
close all
SS = std(STACK,1,3);
pM = bindVec(SS) > graythresh(bindVec(SS))*.2;
imshow(pM,[]);
drawnow
%% make biggest mask for RICH
pM = imdilate(pM,strel('disk',41,0));
imshow(pM,[]);
%% extract SIG from STACK
TOT = max(STACK,[],3) - min(STACK,[],3);
SIG = bsxfun(@minus,STACK,STACK(:,:,1));
SIG = bsxfun(@times,SIG,max(SIG,[],3).^-1);
sz = size(SIG);
SIG = reshape(SIG,[prod(sz(1:2)) sz(3)]);
%{
%% PCA
[SIG C U E L ERR LAM] = PCA_FIT_FULL(double(SIG),3);
%}
%% filter SIG
SIG = imfilter(SIG,ones(1,11)/11,'replicate');
%% reshape SIG
SIG = reshape(SIG,sz);
close all
%{
%% crop for RICH
[J,BOX] = imcrop(STACK(:,:,1),[]);
%}
%% extract wave fronts for MOVIE
close all
%MX = max(dds,[],3);
PER = [.2 .35 .5];
%PER = .35;
VW = zeros([size(STACK,1) size(STACK,2) 3 size(STACK,3)]);
CL = {'r' 'y' 'b'};
disp = 1;
toSave = false;
if toSave
    %v = VideoWriter('/mnt/tetra/nate/arthur2.avi');
    %open(v)
end
SKIP = 1;
for e = 1:SKIP:size(STACK,3)
    
    
   
    if disp
        gI = bindVec(STACK(:,:,e));
        out = bindVec(double(gI)).*pM;

        imshow(out,[]);
        hold on
    end
    for p = 1:numel(PER)
        CON = {};
        dC = contourc(double(SIG(:,:,e).*pM),[PER(p) PER(p)]);
        cnt = 1;
        while ~isempty(dC)
            ptr = dC(:,1);
            CON{cnt} = dC(:,2:1+ptr(2));
            dC(:,1:1+ptr(2)) = [];
            if disp
                plot(CON{cnt}(1,:),CON{cnt}(2,:),CL{p});
                hold on
            end
            cnt = cnt + 1;
           
        end
        STORE{e}{p} = CON;
    end
     drawnow
     hold off
    
    
    %writeVideo(v,getframe(gca));
    
    e
    title(num2str(e))
    drawnow
end
%close(v)
%% view waves from above gather
close all
for e = 1:SKIP:size(STACK,3)
    %gI = bindVec(STACK(:,:,e));
    gI = STACK(:,:,e);
    out = bindVec(double(gI)).*pM;
    imshow(out,[]);
    hold on
    CON = STORE{e};
    for level = 1:numel(CON)
        lCON = CON{level};
        for p = 1:numel(CON)
            plot(lCON{level}(1,:),lCON{level}(2,:),CL{level});
        end
    end
    title(num2str(e))
    drawnow
    hold off
end
%% stack wave fronts for velocity calculation
close all
%MX = max(dds,[],3);
PER = [.2 .35 .5];
PER = .35;
VW = zeros([size(STACK,1) size(STACK,2) 3 size(STACK,3)]);
CL = {'r' 'y' 'b'};
SKIP = 1;
H = {};
tc = 1;
for e = 1:SKIP:size(STACK,3)
    
    
   
    %{
    gI = bindVec(STACK(:,:,e));
    out = bindVec(double(gI)).*FP;
   
    imshow(out,[]);
    hold on
    %}
    for p = 1:numel(PER)
        CON = {};
        dC = contourc(double(SIG(:,:,e).*pM),[PER(p) PER(p)]);
        cnt = 1;
        while ~isempty(dC)
            ptr = dC(:,1);
            CON{cnt} = dC(:,2:1+ptr(2));
            
            
            CON{cnt} = arcLength(CON{cnt}','arcLen');
            
            
            
            TAN = gradient(CON{cnt}');
            L = sum(TAN.*TAN,1).^-.5;
            TAN = bsxfun(@times,TAN,L);
            NOR = [TAN(2,:);-TAN(1,:)];
            NOR = NOR';
            dC(:,1:1+ptr(2)) = [];
            
            
            
            
            H{tc} = [CON{cnt} e*ones(size(CON{cnt},1),1) NOR];
            tc = tc + 1;
           
            
            
            %plot(CON{cnt}(1,:),CON{cnt}(2,:),CL{p});
            %hold on
            cnt = cnt + 1;
            drawnow
        end
        STORE{e}{p} = CON;
    end
    
    
    %writeVideo(v,getframe(gca));
    
    e
    %title(num2str(e))
    %drawnow
end
%close(v)
%% stack fronts for calculation
H = cell2mat(H');
%% 
rm = any(isnan(GGG),2);
GGG(rm,:) = [];
[gS gC gU gE gL gERR gLAM] = PCA_FIT_FULL(GGG,20);
[sweepD] = sweepPCA(gC,gE,gU,std(gC,1,1),1,10);
%% predict velcity
close all
kp = VEL(:,end) < 5;
GGG = GGG(kp,:);
VEL = VEL(kp,:);
[TrainIDX, TestIDX] = crossvalind('HoldOut', size(GGG,1), .2);
 [Xloadings,Yloadings,Xscores,Yscores,  beta,pctVar,mse,stats,Weights] = plsregress(GGG(TrainIDX,:),VEL(TrainIDX,end),20);
 vP = [ones(sum(TestIDX),1) GGG(TestIDX,:)]*beta;
 vPR = [ones(sum(TrainIDX),1) GGG(TrainIDX,:)]*beta;
 plot(vP,VEL(TestIDX,end),'r.')
 hold on
  plot(vPR,VEL(TrainIDX,end),'.');
  %%
  close all
  net = fitnet([5],'trainbr');
  YY = VEL(:,end)';
  net = train(net,gC',YY);
  YYp = net(gC');
  plot(YYp,YY,'.')
%%
for p = 1:size(sweepD,1)
    for step = 1:size(sweepD,2)
        tmpW = sweepD(p,step,:);
        mesh(g1,g2,reshape(tmpW,size(g1,1),size(g1,2)));
        %view([90 90]);
        waitforbuttonpress
    end
end
%%
close all
plot(gC(:,1),gC(:,2),'.')
%% pair wise track fronts
close all
SKIP = 10;
disp = 0;
VEL = [];
stackNUM = 1;
SEL = 0;


[t1,t2] = ndgrid(linspace(0,10,20),linspace(-pi,pi,100));
g1 = t1.*cos(t2);
g2 = t1.*sin(t2);
sampleWAVEfront = 1;
gcnt = 1;
GGG = [];
for e = 1:SKIP:size(STACK,3)
    if disp
        gI = bindVec(STACK(:,:,e));
        out = bindVec(double(gI)).*pM;
        imshow(out,[]);
        hold on
    end
    
    if disp & SEL
        [c r V] = impixel();
    end
   
    if disp
        CON = STORE{e}{1};
        for c = 1:numel(CON)
            plot(CON{c}(:,1),CON{c}(:,2),'r')
        end
    end
    
    
    
    if disp
        CON = STORE{e+SKIP}{1};
        for c = 1:numel(CON)
            plot(CON{c}(:,1),CON{c}(:,2),'b')
        end
    end
    
    
     
    fidx = find(H(:,3)==e);
    
    if disp
        quiver(H(fidx,1),H(fidx,2),-H(fidx,4),-H(fidx,5),.1)
        plot(H(fidx,1),H(fidx,2),'m.')
    end
    
    
    
    fidx1 = find(H(:,3)==e+SKIP);
    if disp
       %  quiver(H(fidx1,1),H(fidx1,2),-H(fidx1,4),-H(fidx1,5),.1)
     %   plot(H(fidx1,1),H(fidx1,2),'c.')
    end
    
    
    if sampleWAVEfront
        if~isempty(fidx) & ~isempty(fidx1)
            tmpI = double(STACK(:,:,e));
            rq = randi(numel(fidx),min(100,numel(fidx)),1);

            
            
            rm = any(isnan(H(fidx(rq),:)) | isinf(H(fidx(rq),:)),2);
            rq(rm) = [];
            
            for f = 1:numel(rq)

                NO = -H(fidx(rq(f)),4:5);
                TN = [-NO(2) NO(1)];
                FRM = [NO',TN'];

                GX = FRM*[g1(:)';g2(:)'];
                GX(1,:) = GX(1,:) + H(fidx(rq(f)),1);
                GX(2,:) = GX(2,:) + H(fidx(rq(f)),2);
                
                %{
                plot(GX(1,:),GX(2,:),'.')
                hold on
                drawnow
                waitforbuttonpress
                %}
                
                GGG(gcnt,:) = ba_interp2(tmpI,GX(1,:),GX(2,:));
                gcnt = gcnt + 1;


            end
        end
    end
    
    tVEL = [];
    if ~isempty(fidx) & ~isempty(fidx1)
        try
            
            
            vi = delaunayTriangulation(H(fidx1,1:2));
            
            %{
            randN = randi(numel(fidx),1000,1);
            fidx1 = fidx(randN);
            %}
            
            fidx = fidx(rq);
            
            nn = [];
            subX = H(fidx,:);
            parfor k = 1:numel(fidx)
                DD = linspace(-20,20,100);
                %searchLN = DD'*H(fidx1(k),4:5);
                %searchLN = bsxfun(@plus,searchLN,H(fidx1(k),1:2));
                searchLN = DD'*subX(k,4:5);
                searchLN = bsxfun(@plus,searchLN,subX(k,1:2));
                %plot(searchLN(:,1),searchLN(:,2),'g.')
                [ntmp,kD] = nearestNeighbor(vi,searchLN);
                [~,sidx] = sort(kD,'ascend');
                [~,midx] = min(abs(DD((sidx(1:10)))));
                nn(k) = ntmp(sidx(midx));
            end
           
            
            
            %nn = nearestNeighbor(vi,H(fidx,1:2));

    
            tVEL = [];
            for pt = 1:numel(nn)
                vect(1,:) = [vi.Points(nn(pt),1) H(fidx(pt),1)];
                vect(2,:) = [vi.Points(nn(pt),2) H(fidx(pt),2)];

                
                tVEL(pt,:) = [H(fidx(pt),:) norm(diff(vect,1,2))];
                if disp
                    if norm(diff(vect,1,2)) < 40
                        plot(vect(1,:),vect(2,:),'y')
                    end
                end
            end
            
            
            VEL = [VEL;tVEL];
            
            if size(VEL,1) ~= size(GGG,1)
                break
            end
        catch ME
            getReport(ME)
        end
    end
    
   
    if disp
        drawnow
    end
    e
end
%% free hand for velocity distribtuions
THRESH = 20;
MAX = 20;
for r = 1:4
    imshow(pM,[]);
    h = imfreehand(gca);
    position = wait(h);
    in = inpolygon(VEL(:,1),VEL(:,2),position(:,1),position(:,2));
    tVEL = VEL(find(in),end);
    tVEL(tVEL > THRESH) = [];
    tVEL(tVEL < 1) = [];
    numel(tVEL)
    [yi{r},xi{r}] = ksdensity(tVEL,linspace(0,MAX,50));
    hold on
end

close all
figure;
CL = {'r' 'g' 'b' ,'k'};
for r = 1:4
    hold on
    yi{r} = yi{r} / sum(yi{r});
    plot(xi{r},yi{r},CL{r})
end
%% area decrease by wave front monster
close all
%MX = max(dds,[],3);
PER = [.2 .35 .5];
PER = .35;
CL = {'r' 'y' 'b'};
SKIP = 1;
FP = pM;
bSTACK = zeros(size(STACK),'single');
for e = 1:SKIP:size(STACK,3)
    
    
    M50 = SIG(:,:,e) < PER;
   
    
    FP = FP.*M50;
    bSTACK(:,:,e) = FP;
    %imshow(FP,[]);
    %drawnow
    %writeVideo(v,getframe(gca));
    
    e
    %title(num2str(e))
    %drawnow
end
%close(v)
%% free hand for velocity distribtuions
close all
THRESH = 20;
MAX = 20;
for r = 1:4
    imshow(pM,[]);
    h = imfreehand(gca);
    position = wait(h);
    areaMSK{r} = poly2mask(position(:,1),position(:,2),size(STACK,1),size(STACK,2));
    imshow(areaMSK{r},[]);
    drawnow
    waitforbuttonpress
    
end

close all

figure;
CL = {'r' 'g' 'b' ,'k'};
areaD = [];
SKIP = 1;
for r = 1:4
    for tm = 1:SKIP:size(bSTACK,3)
        tmp = bSTACK(:,:,tm).*areaMSK{r};
        %imshow(tmp,[]);
        %drawnow
        areaD(r,tm) = sum(tmp(:));
        r
        tm
        %plot(areaD);
        %drawnow
    end
end
%% plot decrease vai wave front and distribution
close all
h1 = figure;
h2 = figure;
h3 = figure;
for r = 1:4
    
    figure(h1);
    sig = areaD(r,:);
    %sig = imfilter(sig,ones(1,51)/51,'replicate');
    velF = gradient(-sig);
    plot(sig,CL{r})
     hold on
     
    figure(h2);
    plot(velF,CL{r})
    hold on
    velF(velF==0) = [];
    
   
    
    %velF(velF > 0) = [];
    
    figure(h3)
    [yi xi] = ksdensity(-velF,linspace(0,200,200));
    plot(xi,yi,CL{r})
    
    hold on
end
%%
close all
UQ = unique(VEL(:,3));
MAX = 20;
xVEL = linspace(0,MAX,100);
for tm = 1:numel(UQ)
    tVEL = VEL(find(VEL(:,3)==UQ(tm)),end);
    tVEL(tVEL > MAX) = [];
    tVEL(tVEL < 1) = [];
    if ~isempty(tVEL)
        YI(tm,:) = ksdensity(tVEL(:),xVEL);
        YI(tm,:)  = YI(tm,:)  / sum(YI(tm,:));
        plot(xVEL,YI(tm,:))
        title(num2str(UQ(tm)));
        drawnow
        waitforbuttonpress
    end
end
%%
close all
mesh(YI);
colormap('jet');
view([0 90]);
%%
close all
%MX = max(dds,[],3);
PER = [.2 .35 .5];
PER = .35;
dPER = .005;
VW = zeros([size(STACK,1) size(STACK,2) 3 size(STACK,3)]);
CL = {'r' 'y' 'b'};
%v = VideoWriter('/mnt/tetra/nate/arthur2.avi');
%open(v)
H = [];
FP = pM;
tc = 1;
for e = 1:2:size(STACK,3)
    
    
    [dF1,dF2] = gradient(STACK(:,:,e));
    
    M50 = SIG(:,:,e) < PER;
   
    
    gI = bindVec(STACK(:,:,e));
    out = bindVec(double(gI)).*FP;
   
    imshow(out,[]);
    hold on
    for p = 1:numel(PER)
        CON = {};
        dC = contourc(double(SIG(:,:,e).*FP),[PER(p) PER(p)]);
        cnt = 1;
        while ~isempty(dC)
            ptr = dC(:,1);
            CON{cnt} = dC(:,2:1+ptr(2));
            
            %ba_interp2(CON{cnt}(
            
            
            %H = [H;[CON{cnt}' e*ones(size(CON{cnt},2),1)]];
            H{tc} = [CON{cnt}' e*ones(size(CON{cnt},2),1)];
            tc = tc + 1;
            dC(:,1:1+ptr(2)) = [];
            
            
            
            plot(CON{cnt}(1,:),CON{cnt}(2,:),CL{p});
            hold on
            cnt = cnt + 1;
        end
        STORE{e} = CON;
    end
    drawnow
    hold off
    
     %FP = FP.*M50;
    
    
    %writeVideo(v,getframe(gca));
    
    %{
    for p = 1:numel(PER)
        MSK = SIG(:,:,e) > (PER(p) - dPER) &  SIG(:,:,e) < (PER(p) + dPER);
        out = flattenMaskOverlay(out,logical(pM.*MSK),1,CL{p});
    end
    %}
    %{
    MSK = SIG(:,:,e) > (PER - dPER) &  SIG(:,:,e) < (PER + dPER);
    MSK = SIG(:,:,e) > (PER - dPER) &  SIG(:,:,e) < (PER + dPER);
    %MSK = SIG(:,:,e) > (.5 - .01);
    MSK2 = SIG(:,:,e) < (.5);
    MSK3 = SIG(:,:,e) > (.5);
    %}
    %gI = bindVec(STACK(:,:,e));
   
    %out = flattenMaskOverlay(out,logical(pM.*MSK2),.1,'r');
    %out = flattenMaskOverlay(out,logical(pM.*MSK3),.1,'b');
    %out= flattenMaskOverlay(out,logical(pM.*MSK2));
    
    %VW(:,:,:,e) = out;
    e
    %imshow(out,[]);
    title(num2str(e))
    drawnow
end
%close(v)
%%
close all
sdF = imfilter(dF,fspecial('disk',21),'replicate');
sdF(isnan(sdF)) = 0;
sdF(isinf(sdF)) = 0;
[v1,v2] = gradient(sdF);
v1(isnan(v1)) = 0;
v1(isinf(v1)) = 0;
v1(abs(v1)>10) = 0;
dt = (v1 + v2);
figure;
imshow(pM.*dt.*(abs(dt)<1),[]);
close all
%%

%%
figure;
imshow(vel < 1,[]);
figure;
dt(abs(dt) > 10) = 0;
imshow((dt).*pM,[]);
%%
close all
imshow(STACK(:,:,1),[]);
for e = 1:15:numel(STORE)
    
    hold on;
    for b = 1:numel(STORE{e})
        plot(STORE{e}{b}(1,:),STORE{e}{b}(2,:),'r')
        
    end
    drawnow
end
%% view
close all
for e = 1:size(VW,4)
    imshow(VW(:,:,:,e),[]);
    title(num2str(e))
    drawnow
end
%%
[c,r,V] = impixel(STACK(:,:,1),[]);
%%
close all
sig = (squeeze(SIG(r,c,:)));
plot(bindVec(sig),'k')
sig = (squeeze(STACK(r,c,:)));
hold on
plot(bindVec(sig))

sig = (squeeze(dds(r,c,:)));
hold on
plot(bindVec(sig))

