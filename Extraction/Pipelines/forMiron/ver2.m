%% read the video from disk
%M = VideoReader('/mnt/tetra/nate/Drug.avi');
M = VideoReader('~/D1-3.hsv');
T = M.NumberofFrames;
H = M.Height;
W = M.Width;
S = zeros(H,W,T,'uint8');
cnt = 1;
%M = VideoReader('/mnt/tetra/nate/Drug.avi');
M = VideoReader('~/D1-3.hsv');
while hasFrame(M)
    tmp = M.readFrame();
    S(:,:,cnt) = tmp(:,:,1);
    cnt = cnt + 1
end
%%
close all
for e = 1:size(S,3)
    imshow(S(:,:,e),[])
    title(num2str(e))
    drawnow
end
%% remove first N frames from stack
S(:,:,1:45) = [];
S(:,:,1700:end) = [];
%% make plate mask
close all
plateMSK = zeros(size(S,1),size(S,2));
cp = hsize(plateMSK);

plateMSK(cp(1),cp(2)) = 1;
plateD = bwdist(plateMSK);
imshow(plateD,[]);
plateMSK = plateD < 500;
imshow(plateMSK,[]);
%% 
S = bsxfun(@times,double(S),plateMSK);
%% get the global stats
N = 10000;
uS = mean(single(S),3);
uSm = min(single(S),[],3);
uSM = max(single(S),[],3);
uRng = uSM - uSm;

sS = std(single(S),1,3);
dS = (diff(single(S),1,3));

udSm = min(single(dS),[],3);
udSM = max(single(dS),[],3);
udRng = udSM - udSm;


dS1 = zeros(size(S));
dS2 = zeros(size(S));
for e = 1:size(S,3)
    [dS1(:,:,e),dS2(:,:,e)] = gradient(single(S(:,:,e)));
    
    e
end

dTOT = (dS1.^2 + dS2.^2).^.5;

dTOTSm = min(single(dTOT),[],3);
dTOTSM = max(single(dTOT),[],3);
dTOTRng = dTOTSM - dTOTSm;

d1Sm = min(single(dS1),[],3);
d1SM = max(single(dS1),[],3);
d1Rng = d1SM - d1Sm;

d2Sm = min(single(dS2),[],3);
d2SM = max(single(dS2),[],3);
d2Rng = d2SM - d2Sm;


sV = std(dS(1:N));
uV = mean(dS(1:N));
%% distribution of window lengths
close all
vec = uRng(:);
vec = vec(randperm(numel(vec)));
vec = vec(1:100000);
options = statset('MaxIter',10);
window1Model = fitgmdist(vec,2,'Options',options);
[c1,~,pc1] = window1Model.cluster(uRng(:));
c1 = reshape(c1,size(uRng));
pc1 = reshape(pc1(:,2),size(uRng));
imshow(c1,[]);
close all
imshow(pc1,[]);
%% distribution of window widths
close all
vec = udRng(:);
vec = vec(randperm(numel(vec)));
vec = vec(1:10000);
ksdensity(vec);
window1Model = fitgmdist(vec,2,'Options',options);
c2 = window1Model.cluster(udRng(:));
c2 = reshape(c2,size(udRng));
imshow(c2,[]);
%% show both
close all
imshow((c1 == 2 & c2 == 2),[]);
%% sample places of interest and record x,y
X = [];
close all

b1 = linspace(0,1,50); % percent in intensity window
b2 = linspace(0,1,20); % percent in velocity window
b3 = linspace(0,1,30); % percet in gradient window
b4 = linspace(0,256,100); % abs location in intensity


statePDF = zeros(numel(b1),numel(b2),numel(b3),numel(b4));



X = [];
for e = 1:10:size(S,3)
    
    
    
N = 500000;
fidx = find((c1 == 2 & c2 == 2));
fidx = fidx(randperm(numel(fidx)));
fidx = fidx(1:N);

    
    q1 = single(S(:,:,e));
    q2 = single(dS(:,:,e));
    q3 = single(dTOT(:,:,e));
    q4 = uRng;
    
    
    s1 = (q1 - uSm).*uRng.^-1;
    s2 = (q2 - udSm).*udRng.^-1;
    s3 = (q3 - dTOTSm).*dTOTRng.^-1;
    s4 = q4;
    
    
    s1 = s1(fidx);
    s2 = s2(fidx);
    s3 = s3(fidx);
    s4 = s4(fidx);
    
    s1 = discretize(s1,b1);
    s2 = discretize(s2,b2);
    s3 = discretize(s3,b3);
    s4 = discretize(s4,b4);
    
    
    idx = sub2ind(size(statePDF),s1,s2,s3,s4);
    
    
    
    idx(isnan(idx)) = [];
    
    [sx,s2,s3,s4] = ind2sub(size(statePDF),idx);
    
    statePDF(idx) = statePDF(idx) + 1;
    
    %X = [X;[sx,sv,sg]];
    
    %vstatePDF = statePDF .* sum(statePDF(:))^-1;
    
    %vstatePDF = imfilter(vstatePDF,fspecial('disk',11),'replicate');
    
    %imshow(interp2(vstatePDF),[]);
    
    %plot(tmpX(:,1),tmpX(:,2),'.');
    %axis([0 1 0 1])
    drawnow
    
    %X = [X;tmpX];
    e
end
%%
B = ones(11,11,11,11);
B = B / sum(B(:));
%statePDFN = convn(statePDF,B,'same');
statePDFN = imfilter(statePDF,B,'replicate');
statePDFN = statePDFN * sum(statePDFN(:))^-1;
%%
statePDFN = statePDF * sum(statePDF(:))^-1;
%%
statePDFN = smooth3(statePDF,'gaussian',[21,21,21],11);
statePDFN = statePDFN * sum(statePDFN(:))^-1;


%{
%%
statePDFN = imfilter(statePDFN,fspecial('disk',11),'replicate');
statePDFN = statePDF * sum(statePDF(:))^-1;
%}
%{
%% find time trends
BX = bsxfun(@minus,S,uSm);
BX = bsxfun(@times,BX,uRng.^-1);

BV = bsxfun(@minus,dS,udSm);
BV = bsxfun(@times,BV,udRng.^-1);

sz1 = size(BX);
sz2 = size(BV);
BX = reshape(BX,[prod(sz1(1:2)) sz1(3)]);
BV = reshape(BV,[prod(sz2(1:2)) sz2(3)]);

kidx1 = ~any(isnan(BX),2);
kidx2 = ~any(isnan(BV),2);

uBX = mean(BX(kidx1,:),1);
uBV = mean(BV(kidx2,:),1);
%%
close all
plot(uBX);
figure;
plot(uBV);
%}
%% make background
warning off
filter_std_set = [5 5];
area_Filter = [20 2400];

perThreshDrop = .9;
diskSize = 11;

% prepare background(s)
BK = double(uS);
SS = imfilter(sS,fspecial('gaussian',[31 31],11),'replicate');

BK = imfilter(uSm,fspecial('gaussian',[31 31],3),'replicate');

%BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(1)),'replicate');
%BK = imopen(BK,strel('disk',diskSize,0));
%BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(2)),'replicate');


% play frame - BK
close all
R = {};

disp = false;
disp2 = true;
clear R
h1 = figure;
h2 = figure;
p1 = [];
p2 = [];
interestMSK = (c1 == 2 & c2 == 2);
X = [];
for e = 1:(size(S,3)-1)
    
    q1 = single(S(:,:,e));
    q2 = single(dS(:,:,e));
    q3 = single(dTOT(:,:,e));
    q4 = uRng;
    
    
    s1 = (q1 - uSm).*uRng.^-1;
    s2 = (q2 - udSm).*udRng.^-1;
    s3 = (q3 - dTOTSm).*dTOTRng.^-1;
    s4 = q4;
    
    s1 = discretize(s1,b1);
    s2 = discretize(s2,b2);
    s3 = discretize(s3,b3);
    s4 = discretize(s4,b4);
    
    
    idx = sub2ind(size(statePDFN),s1,s2,s3,s4);
    
    
    
    idx(isnan(idx)) = 1;
    
    
    Prob = statePDFN(idx);
    Prob = reshape(Prob,size(q1));
    
    %Prob = Prob .* interestMSK;
    
    %Prob(find(interestMSK == 1)) = 1;
    
    mnV = max(Prob(find(interestMSK==1)));
    
    %mini = min(Prob(find(interestMSK==1 & Prob == 0)));
    %Prob(Prob == 0) = mini;
    
    Prob(interestMSK==0) = mnV;
    
    
    Prob = -log(Prob);
    
    
    %Prob(interestMSK==0) = mnV;
    
    
    %Prob = Prob.*interestMSK;
    
    %Prob(isnan(Prob(:))) = 0;
    %Prob(isinf(Prob(:))) = 1;
    
    
    %Prob = bindVec(Prob);
    
    
    
    
    fv = Prob(interestMSK==1);
    
    fv(isinf(fv)) = max(fv(~isinf(fv)));
    
    fv = bindVec(fv);
    th = graythresh(fv);
    th = .08;
    [py,pidx] = ksdensity(fv);
    [py,pidx] = ksdensity(fv,linspace(pidx(1),pidx(end),500));
    bin = imdilate(py,ones(1,21)) == py;
    bidx = find(bin);
    by = py(bidx);
    [~,bb] = sort(by);
    th = 1.5*pidx(bidx(bb(end-1)));
    
    th;
    
    fv = fv > th;
    msk = zeros(size(Prob));
    msk(interestMSK==1) = fv;
    msk = imerode(msk,strel('disk',1,0));
    
  
    
    [msk] = bwareaopenRange(logical(msk),[75 1500]);
    msk = imdilate(msk,strel('disk',1,0));
    R = regionprops(logical(msk),'Centroid','Area','MajorAxisLength','PixelIdxList');
    %R = R([R.MajorAxisLength]<100);
    msk = zeros(size(msk));
    for o = 1:numel(R)
        msk(R(o).PixelIdxList) = 1;
    end
    
    %{
    out = flattenMaskOverlay(double(bindVec(q1)),logical(msk));
    %out = flattenMaskOverlay(out,logical(interestMSK),.8,'b');
    
    imshow(out,[]);
    hold on
    for o = 1:numel(R)
        plot(R(o).Centroid(1),R(o).Centroid(2),'m*');
        %text(R(o).Centroid(1),R(o).Centroid(2),num2str(R(o).MajorAxisLength),'Background','w');
    end
    
    hold off
    
    drawnow
    %}
    
    for o = 1:numel(R)
        X = [X;[R(o).Centroid e]];
    end
    e
    %waitforbuttonpress
end
%{

    
    % mean background subtraction
    tmp = double(S(:,:,e)) - BK;
    tmp = tmp .* (SS+1).^-1;
    
    
    
    
    vtmp = abs((dS(:,:,e) - uV)/sV);
    
    
    
    
    
    
    
    %tmp = imfilter(tmp,fspecial('gaussian',[31 31],11),'replicate');
    % zscore
   
    
    tmp = bindVec(tmp);
    msk = tmp > perThreshDrop*graythresh(tmp);
    vmsk = vtmp > 15;
    
    p1(e) = sum(vmsk(:))/numel(tmp);
    p2(e) = sum(msk(:) | vmsk(:))/numel(tmp);
    
    figure(h2);
    hold on
    plot(p1,'r')
    plot(p2,'b');
    hold off
    
    tmsk = msk | vmsk;
    
    
    
    midx1 = find(tmsk==1);
    midx0 = find(tmsk==0);
    midx0 = midx0(randperm(numel(midx0)));
    midx0 = midx0(1:numel(midx1));
    
    
    if disp
        subplot(1,2,1);
        plot(tmp(midx1),vtmp(midx1),'r.');
        hold on
        plot(tmp(midx0),vtmp(midx0),'b.');
        hold off 
        drawnow
    end

    
    R{e} = regionprops(logical(msk),'Centroid','Area','Perimeter');
    
    
    
    if disp2
        
        figure(h1);
        out = flattenMaskOverlay(bindVec(double(S(:,:,e))),msk,.5,'b');
        out = flattenMaskOverlay(out,vmsk,.5,'r');

        %subplot(1,2,2);

        imshow(out,[]);
        
        
        %imshow(bindVec(tmpv),[]);
        %drawnow
    end
    e
end
%}
%{
        %{
       
        
        
        %tmp = bindVec(tmp);
        msk = tmp > perThreshDrop*graythresh(tmp);
%}
        
        
        
        
        R = regionprops(msk,'Area','Perimeter','Centroid','PixelidxList');
        
        
        %fidx = find([R.Area] > area_Filter(1) & [R.Area] < area_Filter(2));
        %R = R(fidx);
        
        
        msk = zeros(size(msk));
        for o = 1:numel(fidx)
            msk(R(o).PixelIdxList) = 1;
        end

        
        for o = 1:numel(R)
            features{e}(:,o) = [R(o).Centroid';R(o).Area;R(o).Perimeter];
        end
        
        out = flattenMaskOverlay(S(:,:,e),logical(msk));

        imshow(out,[]);
        drawnow
    end
end
%}
%%
D = [X(:,3) X(:,1) X(:,2)];
%% 
D = [D D(:,1)];
D(:,1) = .1*D(:,1);
%% look at some distances
UQ = unique(D(:,end));
for t = 1:(numel(UQ)-1)
    f1 = find(D(:,end) == UQ(t));
    f2 = find(D(:,end) == UQ(t+1));
    s1 = D(f1,1:3);
    s2 = D(f2,1:3);
    
    
    for e = 1:size(s1,1)
        delta = bsxfun(@minus,s2,s1(e,:));
        delta = sum(delta.*delta,2).^.5;
        mdSPACETIME{t}(e) = min(delta);
        
        delta = bsxfun(@minus,s2(:,2:3),s1(e,2:3));
        delta = sum(delta.*delta,2).^.5;
        mdSPACE{t}(e) = min(delta);
        
    end
    
    
    
    MdSPACETIME(t) = max(mdSPACETIME{t});
    MdSPACE(t) = max(mdSPACE{t});
    t
end
%% looking at min distance between objects helps this threshold for connections
R = 20;
A = Radjacency(D(:,1:3)',R);
%% 
A = sparse(size(D,1),size(D,1));
IDX = [];
V = [];
close all
bu = [];
for t = 1:(numel(UQ)-1)
    
    
    f1 = find(D(:,end) == UQ(t));
    f2 = find(D(:,end) == UQ(t+1));
    s1 = D(f1,1:3);
    s2 = D(f2,1:3);
    d = pdist2(s1,s2);
    
    b(t) = min(d(:));
    %d = tril(d);
    zidx = find(d~=0);
    g = d(zidx);
    [g,sidx] = sort(g);
    [idx1,idx2] = ind2sub(size(d),zidx(sidx(1:10)));
    bu(t) = mean(g(1:10));
    %{
    if t == 21
        figure;
        imshow(S(:,:,t),[]);
        hold on
        plot(s1(idx1,2),s1(idx1,3),'r.')
        plot(s2(idx2,2),s2(idx2,3),'go')
        plot(s2(:,2),s2(:,3),'g.')
        pause(1)
        figure;
        imshow(S(:,:,t+1),[]);
        hold on
        plot(s1(idx1,2),s1(idx1,3),'r.')
        plot(s2(idx2,2),s2(idx2,3),'go')
        waitforbuttonpress
        
    end
    %}
    
    d(d > 20) = 0;
    f = cartprod(f1,f2);
    %f = [f(:,2) f(:,1)];
    IDX = [IDX;f];
    V = [V;d(:)];
    t
end 
A = sparse(IDX(:,1),IDX(:,2),V,size(D,1),size(D,1));
%%
for e = 1:numel(UQ)
    fidx = find(D(:,end) == UQ(e));
    for i = 1:numel(fidx)
        for j = 1:numel(fidx)
             A(fidx(i),fidx(j)) = 0;
             A(fidx(j),fidx(i)) = 0;
        end
    end
    e
end
%%
R = sparse(3,3);
R(1,2) = 1;
R(2,3) = 1;
[g1 , g2] = dijkstra(R , 3 , 1);
       
        
%%
close all
fidxStr = find(D(:,end) == 1);
fidxStp = find(D(:,end) == max(UQ));
pathcost = [];
path = {};
close all
for str = 21:numel(fidxStr)
    parfor stp = 1:numel(fidxStp) 
        stp
        [path{str,stp} , pathcost(str,stp)] = dijkstra(A ,fidxStp(stp),fidxStr(str));
       
    end
    
        if ~isempty(path{str,stp})
            for t = 1:10:size(S,3)
                D(path{str,stp}(t),1:4);
                imshow(S(:,:,t),[]);
                hold on
                plot(D(path{str,stp},2),D(path{str,stp},3),'m');
                plot(D(path{str,stp}(t),2),D(path{str,stp}(t),3),'g.');
                hold off
                drawnow
                TAU(t) = D(path{str,stp}(t),4);
                %waitforbuttonpress
                t
            end
            %}
            %close all
            plot(D(:,2),D(:,3),'b.','MarkerSize',1)
            hold on
            plot(D(fidxStr(str),2),D(fidxStr(str),3),'g*','MarkerSize',5)
            plot(D(fidxStp(stp),2),D(fidxStp(stp),3),'ro','MarkerSize',5)
            plot(D(path{str,stp} ,2),D(path{str,stp} ,3),'m','LineWidth',3)
            hold off
            drawnow
            pause(.4)
        else
            imshow(S(:,:,1),[]);hold on
            plot(D(fidxStr(str),2),D(fidxStr(str),3),'g*','MarkerSize',5);
            %figure;
            
            drawnow
            hold off
        end
        
        
        
    end
end
%% extract the X,Y,T
X = [];
for e = 1:numel(R)
    for o = 1:numel(R{e})
        X = [X ; [R{e}(o).Centroid e o]];
    end
end

%% 
[~,BOX] = imcrop(S(:,:,1));
close all
%% make background
bSZ = 21;
filter_std_set = [5 5];
area_Filter = [20 2400];
perThreshDrop = .8;
diskSize = bSZ;

% prepare background(s)
BK = double(uS);
BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(1)),'replicate');
BK = imopen(BK,strel('disk',diskSize,0));
BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(2)),'replicate');

%%
clear J
jBK = double(imcrop(BK,BOX));
for e = 1:size(S,3)
    J(:,:,e) = double(imcrop(S(:,:,e),BOX)) - round(jBK);
    %J(:,:,e) = imfilter(J(:,:,e),fspecial('gaussian',[bSZ bSZ],3),'replicate');
    e
end
%%
clear S
%% try entropy filter
CH = 5;
FR = [];
N = 1200;
MSZ = numel(1:10:(size(J,3)-CH))*N;
%H = zeros(2*256,2*256,MSZ);
yIDX = 200:300;
xIDX = 200:300;
H = zeros(numel(yIDX),numel(xIDX),MSZ);
cnt = 1;
bSZ = 51;
mIDX = 1:(size(J,1)*size(J,2));
mIDX = reshape(mIDX,[size(J,1) size(J,2)]);
mIDX = im2colF(mIDX,[bSZ bSZ],[1 1]);
flag = true;



for e = 1:10:(size(J,3)-CH)
    e
    slice = J(:,:,e:e+CH);
    dSlice = diff(slice,1,3);
    STACK = [];

    %tmp = im2colF(double(slice(:,:,1)),[bSZ bSZ],[1 1]);
    ridx = randperm(size(mIDX,2));

    tH = zeros(2*256,2*256,N);

    for s1 = 1:size(dSlice,3)
        tmp = double(slice(:,:,s1));
        dtmp = double(dSlice(:,:,s1));

        %tmp = im2colF(double(slice(:,:,s)),[bSZ bSZ],[1 1]);
        %dtmp = im2colF(double(dSlice(:,:,s)),[bSZ bSZ],[1 1]);
        %ridx = randperm(size(tmp,2));
        %tmp = tmp(:,ridx);
        %dtmp = dtmp(:,ridx);

        for v = 1:N

            tIDX = mIDX(:,ridx(v));
            %ttH = zeros(256,2*256);
            %HH = tH(:,:,v);
            %ind = sub2ind([size(ttH,1) size(ttH,2)],tmp(:,v)+1,dtmp(:,v)+256);
            ind = sub2ind([size(tH,1) size(tH,2) size(tH,3)],tmp(tIDX)+256,dtmp(tIDX)+256,v*ones(numel(tIDX),1));
            for h = 1:numel(ind)
                %ttH(ind(h)) = ttH(ind(h)) + 1;
                %HH(ind(h)) = HH(ind(h)) + 1;
                tH(ind(h)) = tH(ind(h)) + 1;
            end
            %{
            HH = tH(:,:,v);
            ind = sub2ind([size(H,1) size(H,2)],tmp(:,v)+1,dtmp(:,v)+256);
            for h = 1:numel(ind)
                HH(ind(h)) = HH(ind(h)) + 1;
            end
            tH(:,:,v) = HH;
            %}
            %imshow(HH,[]);
            %drawnow
        end
        %STACK = [STACK;tmp];
    end

    
    
    %{
    for s = 1:N
        imshow(tH(:,:,s),[]);
        drawnow
    end
    %}
    
    H(:,:,cnt:(cnt+N-1)) = tH(yIDX,xIDX,:);
    cnt = cnt + N;
    
    %FR(:,e) = std(STACK,1,1);
end
%%
H(:,:,cnt:end) = [];

%% trim H
xIDX = 200:300;
a = squeeze(mean(sum(H,1),3));
b = squeeze(mean(sum(H,2),3));
plot(a);
%%
close all
yIDX = 200:300;
b = squeeze(mean(sum(H,2),3));
plot(b)
%%
%subH = H(yIDX,xIDX,:);
subH = H;
sz = size(subH);
subH = reshape(subH,[prod(sz(1:2)) sz(3)]);
subH = subH * bSZ.^-2;
%%
[U,E,L] = PCA_FIT_FULL_Tws(subH,3);
%%
C = PCA_REPROJ_T(subH,E,U);
%% sweep
close all
SW = sweepPCA(C',E,U',L.^.5,1,5);
SW = squeeze(SW);
for e = 1:size(SW,1)
    tmp = reshape(SW(e,:),[sz(1:2)]);
    %tmp = reshape(U,[sz(1:2)]);
    %tmp = reshape(SW(e,:)-U',[sz(1:2)]);
    imshow(tmp,[]);
    %waitforbuttonpress
    plot(mean(tmp,2))
    %waitforbuttonpress
end
%%
close all
rE = reshape(bsxfun(@plus,E,U),[sz(1:2) 3]);
imshow(rE(:,:,:),[]);
%%
NG = 3;
options = statset('Display','iter','MaxIter',5);
gmm = fitgmdist(C',NG,'Options',options,'RegularizationValue',.01,'Replicates',10);
%% try entropy filter
CH = 5;
MAS = [];
MAS2 = [];
for e = 1:(size(J,3)-CH)
    slice = J(:,:,e:e+CH);
    dSlice = diff(slice,1,3);
    STACK = [];
    
    tmp = im2colF(double(slice(:,:,1)),[bSZ bSZ],[1 1]);
    tH = zeros(2*256,2*256,size(tmp,2));
    
    for s1 = 1:size(dSlice,3)
        
        %tmp = double(slice(:,:,s));
        %dtmp = double(dSlice(:,:,s));

        
        tmp = im2colF(double(slice(:,:,s1)),[bSZ bSZ],[1 1]);
        dtmp = im2colF(double(dSlice(:,:,s1)),[bSZ bSZ],[1 1]);
        tic
        for v = 1:size(mIDX,2)
            %tIDX = mIDX(:,v);
            %tic
            %ttH = zeros(256,2*256);
            %HH = tH(:,:,v);
            %ind = sub2ind([size(tH,1) size(tH,2)],tmp(:,v)+256,dtmp(:,v)+256);
            %ind = sub2ind([size(tH,1) size(tH,2) size(tH,3)],tmp(mIDX(:,v))+256,dtmp(mIDX(:,v))+256,v*ones(size(tIDX,1),1));
            ind = sub2ind([size(tH,1) size(tH,2) size(tH,3)],tmp(:,v)+256,dtmp(:,v)+256,v*ones(size(tmp,1),1));
            
            %uidx = unique(ind);
            %su = sum(ind==uidx');
            %tH(uidx) = tH(uidx) + su';
            
            for h = 1:numel(ind)
                %ttH(ind(h)) = ttH(ind(h)) + 1;
                %HH(ind(h)) = HH(ind(h)) + 1;
                tH(ind(h)) = tH(ind(h)) + 1;
            end
            
            %tH(:,:,v) = tH(:,:,v) + ttH;
            %tH(:,:,v) = tH(:,:,v) + ttH;
            %tH(:,:,v) = HH;
            %v
            %imshow(HH,[]);
            %drawnow
            %toc
        end
        toc
        s1
        %STACK = [STACK;tmp];
    end
    
    sub_tH = tH(yIDX,xIDX,:);
    sz = size(sub_tH);
    sub_tH = reshape(sub_tH,[prod(sz(1:2)) sz(3)]);
    tC = PCA_REPROJ_T(sub_tH,E,U);
    
    [kidx,nlogl,P] = gmm.cluster(tC');
    
    kidx = col2im(kidx,[bSZ bSZ],[size(slice,1) size(slice,2)]);
    KP = [];
    for k = 1:NG
        KP(:,:,k) = col2im(P(:,k),[bSZ bSZ],[size(slice,1) size(slice,2)]);
        %KP(:,:,k) = bindVec(K(:,:,k));
    end
    
    for k = 1:3
        K(:,:,k) = col2im(tC(k,:),[bSZ bSZ],[size(slice,1) size(slice,2)]);
        %K(:,:,k) = bindVec(K(:,:,k));
    end
    
    MAS(:,:,:,e) = KP;
    MAS2(:,:,:,e) = K;
    e
    cnt
    %FR(:,e) = std(STACK,1,1);
end
%%
close all
h1 = figure;
h2 = figure;
for e = 1:size(MAS,4)
   
    [~,midx] = max(MAS(:,:,:,e),[],3);
    msk = midx == 3;
   
    figure(h2)
    tmp = J(:,:,e);
    
    for r = 1:4
        tmp(1:25,:) = [];
        tmp = imrotate(tmp,90);
    end
    imshow(tmp,[])
    
    out = flattenMaskOverlay(bindVec(tmp),msk);
    figure(h1)
    imshow(out,[]);
    drawnow
    %waitforbuttonpress
end
%% get the mean
uS = mean(single(S),3);
%% match features
for e = 1:numel(features)
    f1 = features{e};
    f2 = features{e+1};
    D = zeros(size(f1,2),size(f2,2));
    for i = 1:size(f1,2)
        for j = 1:size(f2,2)
            delta = f1(:,i) - f2(:,j);
            D(i,j) = sum(delta(1:2).*delta(1:2),1).^.5;
        end
    end
    
    for i = 1:size(f1,2)
        [r,sidx] = sort(D(i,:));
        fgamma{e}(i,:) = sidx(1:3);
    end
    e
end
%% plot gamma
P = [];
for e = 1:size(S,3)
    if e == 1
        for pt = 1:size(features{e},2)
            P(pt,:,e) = features{e}(1:2,pt)';
        end
        
    else
        
        for pt = 1:size(P,1)
            P(pt,:,e) = features{e}(1:2,fgamma{e-1}(pt,1))';
        end
    end
    
    
    imshow(S(:,:,e),[]);
    hold on
    for pt = 1:size(P,1)
        trace = squeeze(P(pt,:,:));
        %plot(trace(1,:),trace(2,:),'r')
    end
    
    e
    
end
%% remove the mean frame from the stack
sS = std(single(S),1,3);
Zthreshold = 2;
for e = 1:size(S,3)
    MSK = (single(S(:,:,e)) - uS) > Zthreshold*sS;
    %S1(:,:,e) = double(MSK).*double(S(:,:,e)) + double(~MSK).*uS;
    S1(:,:,e) = double(MSK).*double(S(:,:,e));
    %imshow(S1(:,:,e),[]);
    %title(num2str(e))
    %drawnow
    %imshow(S1(:,:,e),[]);
    %drawnow
    e
    size(S,3)
end
%% make dt and ds
%S1 = bsxfun(@minus,single(S),uS);
DT = diff(S1,1,3);
for e = 1:size(S1,3)
    [d1,d2] = gradient(S1(:,:,e));
    DS(:,:,e) = (d1.^2 + d2.^2).^.5;
    e
end
%%
%{
%%
for e = 1:10%size(S,3)
    tmp = imfilter(S(:,:,e),fspecial('gaussian',[31 31],7),'replicate');
    sz = size(tmp);
    tmp = imresize(tmp,.25);
    tmp = imopen(tmp,strel('disk',7,0));
    tmp = imresize(tmp,sz);
    tmp = imfilter(tmp,fspecial('gaussian',[31 31],7),'replicate');
    tmp = S(:,:,e) - tmp;
    
    pS(:,:,e) = tmp;
    imshow(pS(:,:,e),[]);
    drawnow
end
%%
tmp = double(mean(S,3));
tmp = imfilter(tmp,fspecial('gaussian',[31 31],7),'replicate');
sz = size(tmp);
tmp = imresize(tmp,.25);
tmp = imopen(tmp,strel('disk',7,0));
tmp = imresize(tmp,sz);
tmp = imfilter(tmp,fspecial('gaussian',[31 31],7),'replicate');
BK = tmp;
for e = 1:10%size(S,3)
    pS(:,:,e) = double(S(:,:,e)) - BK;
    imshow(pS(:,:,e),[]);
    drawnow
end
%%
sD = diff(S,1,3);
%%
S = S(:,:,1:end-1);
%}
%% gather the feature space over some of the image patches
% image patches are drawn from S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patchSKUIP = 4;
W = [];
for e = 1:patchSKUIP:size(S1,3)
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = S1(:,:,e);
    H0 = im2colF(double(tmp),[51 51],[11 11]);
    idx = randperm(size(H,2));
    H0 = H0(:,idx);
    H0 = H0(:,1:5000);
    H0 = sort(H0,1);
    H0 = H0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = DT(:,:,e);
    H1 = im2colF(double(tmp),[51 51],[11 11]);
    %idx = randperm(size(H,2));
    H1 = H1(:,idx);
    H1 = H1(:,1:5000);
    H1 = sort(H1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = DS(:,:,e);
    H2 = im2colF(double(tmp),[51 51],[11 11]);
    %idx = randperm(size(H,2));
    H2 = H2(:,idx);
    H2 = H2(:,1:5000);
    H2 = sort(H2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cat the feature vectors together
    W = [W [H0;H1;H2]];
    
    %W = [W H];
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA on the data each seperate feature space - [X,dS,dT];
nD = 3;
%nD = 1;
CHUNK = size(W,1)/nD;
IDX1 = 1:CHUNK;
IDX2 = (CHUNK+1):(2*CHUNK);
IDX3 = (2*CHUNK+1):(3*CHUNK);
[U1,E1,L1] = PCA_FIT_FULL_Tws(W(IDX1,:),3);
[U2,E2,L2] = PCA_FIT_FULL_Tws(W(IDX2,:),3);
[U3,E3,L3] = PCA_FIT_FULL_Tws(W(IDX3,:),3);
L1 = L1(1).^.5;
L2 = L2(1).^.5;
L3 = L3(1).^.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% back project the data
Sc1 = PCA_REPROJ_T(W(IDX1,:),E1,U1);
Sc2 = PCA_REPROJ_T(W(IDX2,:),E2,U2);
Sc3 = PCA_REPROJ_T(W(IDX3,:),E3,U3);
Sc1 = bsxfun(@times,Sc1,L1^-1);
Sc2 = bsxfun(@times,Sc2,L2^-1);
Sc3 = bsxfun(@times,Sc3,L3^-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% direct sum the data together and perform PCA across spaces
ScT = [Sc1;Sc2;Sc3];
[UT,ET,LT] = PCA_FIT_FULL_Tws(ScT,size(ScT,1));
SccT = PCA_REPROJ_T(ScT,ET,UT);
close all
plot(cumsum(LT)/sum(LT));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the PCA of the raw data - space together
sel = nchoosek(1:6,2);
for e = 1:size(sel,1)
    plot(SccT(sel(e,1),:),SccT(sel(e,2),:),'.');
    drawnow
    waitforbuttonpress
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the PCA of the raw data - each space seperate
close all
figure
plot(Sc1(1,:),Sc1(2,:),'.');
figure
plot(Sc1(1,:),Sc1(3,:),'.');
figure
plot(Sc2(1,:),Sc2(2,:),'.');
figure
plot(Sc3(1,:),Sc3(2,:),'.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% "straighten" the data
funcS{1} = @(value)min(value);
funcS{2} = @(value)max(value);
funcS{3} = @(value)mean(value);
funcS{4} = @(value)std(value);
funcS{5} = @(value).5*(max(value) - min(value)) + min(value);

% create combine function of the functions run on the data
combineFunction = @(data,resultY)(data - resultY(5,:));
combineFunction = @(data,resultY)(data - resultY(1,:));

%%%%%%%%%%%%%%%%
% make PC2 a function of PC1
% Data = X space
NSP = 200;
[funcY1] = straightMiron(Sc1,[1 2],funcS,NSP);
nSc1 = applyStraightenMiron(Sc1,funcY1,combineFunction);
close all
plot(Sc1(1,:),Sc1(2,:),'.')
hold on
plot(nSc1(1,:),nSc1(2,:),'r.')

%%%%%%%%%%%%%%%%
% make PC2 a function of PC1
% Data = dS space
NSP = 200;
[funcY2] = straightMiron(Sc2,[1 2],funcS,NSP);
nSc2 = applyStraightenMiron(Sc2,funcY2,combineFunction);
figure
plot(Sc2(1,:),Sc2(2,:),'.')
hold on
plot(nSc2(1,:),nSc2(2,:),'r.')


%%%%%%%%%%%%%%%%
% make PC2 a function of PC1
% Data = dT space
NSP = 200;
[funcY3] = straightMiron(Sc3,[1 2],funcS,NSP);
nSc3 = applyStraightenMiron(Sc3,funcY3,combineFunction);
figure
plot(Sc3(1,:),Sc3(2,:),'.')
hold on
plot(nSc3(1,:),nSc3(2,:),'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% investigate the relationship between PC1 and PC2
K = convhulln(Sc1');
%{
sweepNP = 100;
sweepDIM = 1;
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve(1:2,:) = [sweepDomain;funcY1(1).func(sweepDomain)];
%}
%{
close all
plot(dataToSweep(1,:),dataToSweep(2,:),'.')
hold on
plot(sweepCurve(1,:),sweepCurve(2,:),'r')
waitforbuttonpress
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
close all
plot(sweepF)
%}

close all
figure
plot(Sc1(1,:),Sc1(3,:),'.')
figure
plot(Sc1(1,:),Sc1(2,:),'.')
hold on
plot(Sc1(1,:),nSc1(2,:),'.')
figure
plot(Sc1(2,:),Sc1(3,:),'.')
hold on
plot(nSc1(2,:),Sc1(3,:),'.')
figure;
plot3(Sc1(1,:),Sc1(2,:),Sc1(3,:),'.')
hold on;
trisurf(K,Sc1(1,:)',Sc1(2,:)',Sc1(3,:)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sweep a select line in PC space to see what it does to the eigen function
%% first sweep the orginal - PC1
dataToSweep = Sc1;
sweepDIM = 1;
sweepNP = 10;
sweepRange = std(Sc1,1,2);
sweepDomain = linspace(-sweepRange(sweepDIM),sweepRange(sweepDIM),sweepNP);
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve = zeros(size(dataToSweep,1),sweepNP);
sweepCurve(sweepDIM,:) = L1*sweepDomain;
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
plot(sweepF)
waitforbuttonpress
%% first sweep the orginal - PC2
dataToSweep = Sc1;
sweepDIM = 2;
sweepNP = 10;
sweepRange = std(Sc1,1,2);
sweepDomain = linspace(-sweepRange(sweepDIM),sweepRange(sweepDIM),sweepNP);
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve = zeros(size(dataToSweep,1),sweepNP);
sweepCurve(sweepDIM,:) = L1*sweepDomain;
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
plot(sweepF)
waitforbuttonpress
%% look at the min curve
sweepDIM = 1;
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve(1:2,:) = [sweepDomain;funcY1(1).func(sweepDomain)];
close all
plot(dataToSweep(1,:),dataToSweep(2,:),'.')
hold on
plot(sweepCurve(1,:),sweepCurve(2,:),'r')
waitforbuttonpress
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
close all
plot(sweepF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the crosses between [X,dS,dT]
close all
plot(nSc1(1,:),nSc3(1,:),'.')
figure
plot(nSc3(1,:),nSc1(1,:),'.')
figure;
plot(nSc1(2,:),nSc3(2,:),'.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% twist the crosses - make the PC's of an element of [X,dS,dT] a function of the other
close all
NSP = 500;
dataX1 = cat(1,nSc3(1,:),nSc1(1,:));
[funcY3_cross1] = straightMiron(dataX1,[1 2],funcS,NSP);
ndataX1 = applyStraightenMiron(dataX1,funcY3_cross1,combineFunction);
[funcY3_cross2] = straightMiron(ndataX1,[2 1],funcS,NSP);
ndataX2 = applyStraightenMiron(ndataX1,funcY3_cross2,combineFunction);
figure
plot(dataX1(1,:),dataX1(2,:),'.')
hold on
figure
plot(ndataX1(1,:),ndataX1(2,:),'r.')
figure;
plot(ndataX2(1,:),ndataX2(2,:),'g.')
[Z] = myPDF2(ndataX1,[100 100]);
figure;
lZ = -log(Z);
lZ(isinf(lZ(:))) = 0;
imshow(lZ,[]);
%{
%%%%%%%%%%%%%%%%%%%%%%
% NOTE THAT GMM @ THIS POINT SEEM TO BE A BIT EARLY - 
% RATHER NEEDING TO EXPLORE - AS OF SEP 25,  2019
%%%%%%%%%%%%%%%%%%%%%%
%% fit GMM model(s) to the data
close all
NG = 2;
dataToFit = ndataX1;
kidx = kmeans(dataToFit',NG);

options = statset('Display','iter','MaxIter',1000);
gm = fitgmdist(dataToFit',2,'Options',options,'Replicates',2,'Start','plus');
gmPDF = @(x,y)pdf(gm,[x y]);
figure;
CL = {'r.','g.','b.'};
hold on
for k = 1:NG
    %scatter(dataToFit(1,kidx==k),dataToFit(2,kidx==k),10,CL{k})
end
kidx2 = gm.cluster(dataToFit');
hold on
for k = 1:NG
    scatter(dataToFit(1,kidx2==k),dataToFit(2,kidx2==k),10,CL{k})
end
[n1,n2] = ndgrid(linspace(min(dataToFit(1,:)),max(dataToFit(1,:)),100),linspace(min(dataToFit(2,:)),max(dataToFit(2,:)),100));
E = gmPDF(n1(:),n2(:));
E = reshape(E,size(n1));
hold on
contour(gca,n1,n2,E,10);
hold on
idxL = kidx2;
BKIDX = mode(idxL);
%}
%% second level PCA
subIDX = idxL~=BKIDX;
CHUNK = size(W,1)/3;
IDX1 = 1:CHUNK;
IDX2 = (CHUNK+1):(2*CHUNK);
IDX3 = (2*CHUNK+1):(3*CHUNK);
[U1_2,E1_2,L1_2] = PCA_FIT_FULL_Tws(W(IDX1,subIDX),3);
[U2_2,E2_2,L2_2] = PCA_FIT_FULL_Tws(W(IDX2,subIDX),3);
[U3_2,E3_2,L3_2] = PCA_FIT_FULL_Tws(W(IDX3,subIDX),3);
L1_2 = L1_2(1);
L2_2 = L2_2(1);
L3_2 = L3_2(1);
%% back project the data
Sc1_2 = PCA_REPROJ_T(W(IDX1,subIDX),E1_2,U1_2);
Sc2_2 = PCA_REPROJ_T(W(IDX2,subIDX),E2_2,U2_2);
Sc3_2 = PCA_REPROJ_T(W(IDX3,subIDX),E3_2,U3_2);
Sc1_2 = bsxfun(@times,Sc1_2,L1_2^-1);
Sc2_2 = bsxfun(@times,Sc2_2,L2_2^-1);
Sc3_2 = bsxfun(@times,Sc3_2,L3_2^-1);
%% look at the PCA of the raw data
close all
figure
plot(Sc1_2(1,:),Sc1_2(2,:),'.');
figure
plot(Sc2_2(1,:),Sc2_2(2,:),'.');
figure
plot(Sc3_2(1,:),Sc3_2(2,:),'.');
%% "straighten" the data after section out the data
funcS{1} = @(value)min(value);
funcS{2} = @(value)max(value);
funcS{3} = @(value)mean(value);
funcS{4} = @(value)std(value);
combineFunction = @(data,resultY)(data - resultY(1,:));

NSP = 200;
[funcY1_2] = straightMiron(Sc1_2,[1 2],funcS,NSP);
nSc1_2 = applyStraightenMiron(Sc1_2,funcY1_2,combineFunction);
close all
plot(Sc1_2(1,:),Sc1_2(2,:),'.')
hold on
plot(nSc1_2(1,:),nSc1_2(2,:),'r.')


NSP = 200;
[funcY2_2] = straightMiron(Sc2_2,[1 2],funcS,NSP);
nSc2_2 = applyStraightenMiron(Sc2_2,funcY2_2,combineFunction);
figure
plot(Sc2_2(1,:),Sc2_2(2,:),'.')
hold on
plot(nSc2_2(1,:),nSc2_2(2,:),'r.')



NSP = 200;
[funcY3_2] = straightMiron(Sc3_2,[1 2],funcS,NSP);
nSc3_2 = applyStraightenMiron(Sc3_2,funcY3_2,combineFunction);
figure
plot(Sc3_2(1,:),Sc3_2(2,:),'.')
hold on
plot(nSc3_2(1,:),nSc3_2(2,:),'r.')

%% look at the crosses
close all
plot(nSc1_2(1,:),nSc3_2(1,:),'.')
figure
plot(nSc3_2(1,:),nSc1_2(1,:),'.')
figure;
plot(nSc1_2(2,:),nSc3_2(2,:),'.')
%% twist the crosses
close all
NSP = 500;
dataX1 = cat(1,nSc3_2(1,:),nSc1_2(1,:));
[funcY3_cross1] = straightMiron(dataX1,[1 2],funcS,NSP);
ndataX1 = applyStraightenMiron(dataX1,funcY3_cross1,combineFunction);
[funcY3_cross2] = straightMiron(ndataX1,[2 1],funcS,NSP);
ndataX2 = applyStraightenMiron(ndataX1,funcY3_cross2,combineFunction);
figure
plot(dataX1(1,:),dataX1(2,:),'.')
hold on
figure
plot(ndataX1(1,:),ndataX1(2,:),'r.')
figure;
plot(ndataX2(1,:),ndataX2(2,:),'g.')
[Z] = myPDF2(ndataX1,[100 100]);
figure;
lZ = -log(Z);
lZ(isinf(lZ(:))) = 0;
imshow(lZ,[]);

%% fit GMM model(s) to the data
close all
NG = 2;
dataToFit = ndataX1;
kidx = kmeans(dataToFit',NG);

options = statset('Display','iter','MaxIter',1000);
gm = fitgmdist(dataToFit',2,'Options',options,'Replicates',2,'Start','plus');
gmPDF = @(x,y)pdf(gm,[x y]);
figure;
CL = {'r.','g.','b.'};
hold on
for k = 1:NG
    %scatter(dataToFit(1,kidx==k),dataToFit(2,kidx==k),10,CL{k})
end
kidx2 = gm.cluster(dataToFit');
hold on
for k = 1:NG
    scatter(dataToFit(1,kidx2==k),dataToFit(2,kidx2==k),10,CL{k})
end
[n1,n2] = ndgrid(linspace(min(dataToFit(1,:)),max(dataToFit(1,:)),100),linspace(min(dataToFit(2,:)),max(dataToFit(2,:)),100));
E = gmPDF(n1(:),n2(:));
E = reshape(E,size(n1));
hold on
contour(gca,n1,n2,E,10);
hold on
idxL = kidx2;
BKIDX = mode(idxL);

%% look at scatter plot of sc1
close all
plot(Sc1_2(1,:),Sc1_2(2,:),'.')
figure;
close all
plot(Sc1_2(2,:),Sc2_2(2,:),'.')
%%
%% second level cluster
close all
grps2 = 4;
GMModel_2 = fitgmdist([Sc1_2;Sc2_2]',grps2);
[idxL_2] = cluster(GMModel_2,[Sc1_2;Sc2_2]');
CL = {'r.','g.','b.','k.'};
h1 = figure;
h2 = figure;
for u = 1:grps2
    sidx = idxL_2 == u;
    figure(h1);
    plot(Sc1_2(1,sidx),Sc2_2(1,sidx),CL{u});
    hold on
    figure(h2)
    ksdensity(Sc1_2(1,sidx))
    hold on
end
%%
UQ = unique(idxL);
selIDX = idxL ~= BKIDX;
GMModel2 = fitgmdist([Sc1(:,selIDX);Sc2(:,selIDX)]',3);
%%
plot(sum(W(IDX1,:),1),Sc1(1,:),'.');
plot(Sc1(1,:),Sc1(2,:),'.');
%%
ds = 5;
cnt = 1;
vec1 = [];
STORE = [];
TARN = [];
NN = [];
cnt = 1;
numIMG = size(S1,3)
numIMG = 100;
%numIMG = 10;
R = 15;
N = 25;
MSZ = 500;
MASK = zeros(size(S1,1),size(S1,2));
MASK(end/2,end/2) = 1;
MASK = bwdist(MASK) < MSZ;
for e = 1:1:numIMG
    tic
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = S1(:,:,e);
    reOLD = size(tmp);
    
    NN(e) = randi([15,30],1,1);
    tmp = generateImageP(MASK,R,NN(e));
    imshow(tmp,[]);
    drawnow
    H1 = im2colF(double(tmp),[51 51],[ds ds]);
    
    
    
    
    %H1 = im2colF(double(tmp),[51 51],[(51-1) (51-1)]);
    H1 = sort(H1,1);
    
    % generte targets for detectors
    TARN(cnt,:) = sum(H1,1)*(pi*R)^-2;
    
    
    
    Sc1 = PCA_REPROJ_T(H1,E1,U1);
    Sc1 = bsxfun(@times,Sc1,L1.^-1);
    
    %Sc1_2 = PCA_REPROJ_T(H1,E1_2,U1_2);
    %Sc1_2 = bsxfun(@times,Sc1_2,L1_2.^-1);
    
    STORE(:,:,cnt) = Sc1;
    cnt = cnt + 1;
    
    %[Count,img] = countOP(Sc1,x,25);
    % out = 
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = D1(:,:,e);
    H2 = im2colF(double(tmp),[51 51],[ds ds]);
    H2 = sort(H2,1);
    
    Sc2 = PCA_REPROJ_T(H2,E2,U2);
    Sc2 = bsxfun(@times,Sc2,L2.^-1);
    
    Sc2_2 = PCA_REPROJ_T(H2,E2_2,U2_2);
    Sc2_2 = bsxfun(@times,Sc2_2,L2_2.^-1);
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diff data
    tmp = DS(:,:,e);
    reOLD = size(tmp);
    H3 = im2colF(double(tmp),[51 51],[ds ds]);
    H3 = sort(H3,1);
    
    Sc3 = PCA_REPROJ_T(H3,E3,U3);
    Sc3 = bsxfun(@times,Sc3,L3.^-1);
    Sc3_2 = PCA_REPROJ_T(H3,E3_2,U3_2);
    Sc3_2 = bsxfun(@times,Sc3_2,L3_2.^-1);
    
    vec2(cnt,:) = mean(Sc3,2)';
    cnt = cnt + 1;
    
    
    
    %[idx1,~,P1] = cluster(GMModel,[Sc1;Sc2;Sc3]');
    [idx1,~,P1] = cluster(GMModel,[Sc1]');
    [idx2,~,P2] = cluster(GMModel_2,[Sc1_2;Sc2_2]');
   
    for k = 1:size(P1,2)
        IMG1(:,:,k,e) = col2im(P1(:,k),[51 51],size(tmp));
        IMG_LAB1(:,:,e) = col2im(idx1,[51 51],size(tmp));
    end
    
    
    for k = 1:size(P2,2)
        IMG2(:,:,k,e) = col2im(P2(:,k),[51 51],size(tmp));
        IMG_LAB2(:,:,e) = col2im(idx2,[51 51],size(tmp));
    end
    
    vec1(cnt,:) = mean(Sc1,2)';
    %}
    
    %{
    % back apply beta
    PM = beta'*[Sc1;Sc3];
    PM = col2im(PM,[51 51],size(tmp));
    %}
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = D1(:,:,e);
    H = im2colF(double(tmp),[51 51],[ds ds]);
    H = sort(H,1);
    Sc2 = PCA_REPROJ_T(H,E2,U2);
    Sc2 = bsxfun(@times,Sc2,L2.^-.5);
   
    
    WSc = [Sc1;Sc2];
    [idx,nlogl,P] = cluster(GMModel,WSc');
    
    [idx2,nlogl2,P2] = cluster(GMModel2,WSc');
    
    
    
    LAB1 = col2im(idx,[51 51],size(tmp));
    for k = 1:size(P,2)
        tmp1(:,:,k) = col2im(P(:,k),[51 51],reOLD);
    end
    
    
    LAB2 = col2im(idx2,[51 51],size(tmp));
    for k = 1:size(P2,2)
        tmp2(:,:,k) = col2im(P2(:,k),[51 51],reOLD);
    end
    
    
    
    IMG(:,:,:,e) = bsxfun(@times,tmp2,LAB1~=BKIDX);
    
    
    
    imshow(IMG(:,:,:,e),[]);
    drawnow
    toc
    %}
    toc
end
%%
func = @(vec)countOP(STORE,vec,25);
options = optimoptions('fmincon','Display','iter');
x0 = 2*(rand(1,5)-.5);
x0(end) = 10000/2;
x = fmincon(func,x0,[],[],[],[],[-1;-1;-1;-100;10^3],[1;1;1;100;10^4],[],options);
%%
options = optimset('Display','iter');
x = fminsearch(func,2*(ones(1,5)-.5),options);
%%
func = @(vec)countOP(STORE,vec,25);
options = optimoptions('particleswarm','Display','iter');
x = particleswarm(func,5,[-1;-1;-1;-100;10^3],[1;1;1;100;10^4],options);
%%
func = @(vec)countOP(STORE,vec,25);
options = optimoptions('particleswarm','Display','iter');
x = particleswarm(func,5,[-1;-1;-1;-100;10^3],[1;1;1;100;10^4],options);
%%
sz = [3 3];
func = @(vec)countOP3(STORE,vec,25,sz);
options = optimoptions('particleswarm','Display','iter','MaxIter',200);
xo = [rand(1,prod(sz)) 0 0 0  10 10 10];
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
x = particleswarm(func,numel(xo),LB,UB,options);
%%
global BD
BD = [];
sz = [3 3];
func = @(vec)countOP4(STORE,vec,NN,TARN,[],sz);
OF = @(A,B)outputPLT(func,A,B);
options = optimoptions('particleswarm','Display','iter','MaxIter',200,'UseParallel',true,'PlotFcn',OF);
UU = mean(mean(STORE,2),3);
SU = std(mean(STORE,2),1,3);
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
LB = [-ones(prod(sz),1);UU-10*SU;SU*.0001];
UB = [ones(prod(sz),1);UU+10*SU;SU/2];
xo = [rand(1,prod(sz)) 0*ones(1,sz(1)) 1*ones(1,sz(1))];
x = particleswarm(func,numel(xo),LB,UB,options);
%%
global BD
BD = [];
sz = [3 3];
func = @(vec)countOP4(STORE,vec,NN,TARN,[],sz);
OF = @(A,B,C)outputPLT(func,A,B,C);
options = optimoptions('fmincon','Display','iter','PlotFcn',OF,'UseParallel',true);
xo = [rand(1,prod(sz)) 0*ones(1,sz(1)) 1*ones(1,sz(1))];
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
x = fmincon(func,xo,[],[],[],[],LB,UB,[],options);
%% try count
tmp = S1(:,:,e);
for e = 1:10
    %{
    tmpN = randi([15,30],1,1);
    tmp = generateImageP(MASK,R,tmpN);
    H1 = im2colF(double(tmp),[51 51],[ds ds]);
    H1 = sort(H1,1);
    Sc1 = PCA_REPROJ_T(H1,E1,U1);
    Sc1 = bsxfun(@times,Sc1,L1.^-1);
    %}
    Sc1 = STORE(:,:,e);
    tmpN = NN(e);
    [~,~,tmpCNT] = countOP4(Sc1,x,25,[],[],sz);
    tmpCNT
    tmpN
end
%% next with real data - not yet
global BD
BD = [];
sz = [3 3];
func = @(vec)countOP4(STORE,vec,25,TARN,[],sz);
OF = @(A,B,C)outputPLT(func,A,B,C);
options = optimoptions('fmincon','Display','iter','PlotFcn',OF);
xo = [rand(1,prod(sz)) 0*ones(1,sz(1)) 1*ones(1,sz(1))];
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
x = fmincon(func,xo,[],[],[],[],LB,UB,[],options);
%%
[CNT,IMG] = countOP3(SCORE,x,25,sz);
for k = 1:size(IMG,1)
    IMGV(:,:,k) = col2im(IMG(k,:),[51 51],size(tmp));
end
%%
func = @(vec)countOP2(STORE,vec,25);
T = reshape(STORE,[size(STORE,1) size(STORE,2)*size(STORE,3)]);
kidx = kmeans(T',2);
lda = myLDA(T',kidx);
scores = lda'*T;
for u = 1:2
    U(u) = mean(scores(kidx==u));
    STD(u) = std(scores(kidx==u));
end
vec = [[U STD] lda'];
options = optimoptions('particleswarm','Display','iter');
x = fminsearch(func,vec);
%% beta for counts on vec1 only
beta = vec1\ones(size(vec1,1),1);
%% beta for counts on vec2 and vec1
VEC = [vec1 vec2];
beta = VEC\ones(size(VEC,1),1);

%%
net = feedforwardnet(1);
view(net)

%%
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL_T(W,3);
GMModel = fitgmdist(wC',3);
%%
wSZ = [51 51];
for e = 1%1:1:size(S,3)
    tic
    %tmp = S(:,:,e);
    [tmp] = imcrop(S(:,:,e),BOX);
    %mS(:,:,e) = imresize(S(:,:,e),1);
    H = im2colF(double(tmp)/255,wSZ,[1 1]);
    H1 = sort(H,1);
    
    [tmp] = imcrop(sD(:,:,e),BOX);
    %mS(:,:,e) = imresize(S(:,:,e),1);
    H = im2colF(double(tmp)/255,wSZ,[1 1]);
    H = sort(H,1);
    
    
    H = PCA_REPROJ_T([H;H1],wE,wU);
    [IDX,nlogl,P] = GMModel.cluster(H');
    IDX = col2im(IDX,wSZ,size(tmp));
    PK = col2im(P(:,3),wSZ,size(tmp));
end
%%
close all
for f = 1:size(S,3)
    imshow(S(:,:,f),[]);
    title(num2str(f))
    drawnow
end
%%
S = double(S)/255;
%%
close all
MASK = S(:,:,1) > graythresh(S(:,:,1))*.8;
imshow(MASK,[]);
drawnow
figure;
imshow(S(:,:,1),[]);
fidx = find(MASK);
%%
close all
imshow(mean(S,3),[]);
%% 
S0 = S(:,:,45:end);
uS = mean(S0,3);
uS = imfilter(uS,fspecial('gaussian',[51 51],9));
S0 = bsxfun(@minus,S0,uS);
imshow(uS,[]);
%%
close all
[n1,n2] = ndgrid(100:20:1000,150:20:1100);
pointList = [n1(:) n2(:)];
[pointListL pointListE SE SL] = wholeTrack_mod0(S0,pointList,1,[]);
%%
close all
for e = 1:size(S0,3)
    imshow(S0(:,:,e),[]);
    hold on
    plot(pointListL(:,2,e),pointListL(:,1,e),'r.')
    hold off
    drawnow
end
%%
pointListEM = pointListE;
pointListEM(:,:,1) = 0;
%%
clear X
[X(:,2) X(:,1) v] = impixel(S0(:,:,1),[]);
for t = 1:size(pointListE,3)
    F1{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,1,t),size(n1)));
    F2{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,2,t),size(n1)));
end
close all
clear X
MASK = S(:,:,1) > graythresh(S(:,:,45));
[X(:,1) X(:,2)] = find(MASK);
X = [n1(:) n2(:)];
for t = 1:size(pointListE,3)
    imshow(S0(:,:,t),[]);
    hold on
    plot(X(:,2,t),X(:,1,t),'r.')
    hold off
    drawnow
    d1 = F1{t}(X(:,1),X(:,2));
    d2 = F2{t}(X(:,1),X(:,2));
    X(:,1,t+1) = X(:,1,t) + d1;
    X(:,2,t+1) = X(:,2,t) + d2;
end
%%
close all
SKIP = 15;
[m1,m2] = ndgrid(100:SKIP:1000,150:SKIP:1100);
X = [m1(:) m2(:)];
d = 0;
for t = 1:size(pointListE,3)
    F1{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,1,t),size(n1)),'spline');
    F2{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,2,t),size(n1)),'spline');
end

SKIP = 1;
[xm1,xm2] = ndgrid(100:SKIP:1000,150:SKIP:1100);
XM = [xm1(:) xm2(:)];
d = 0;
for t = 1:size(pointListE,3)
    F1{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,1,t),size(n1)),'spline');
    F2{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,2,t),size(n1)),'spline');
end


PM = zeros(size(xm1));
oX = X;
oXM = XM;
V = [];
VELT = 1.5;



h1 = figure;
h2 = figure;
dataI = zeros(size(xm2,1),size(xm2,2),size(pointListE,3),3);
for t = 1:size(pointListE,3)
    Y = ba_interp2(S(:,:,t+44),xm2,xm1);
    YW = S(:,:,t+45);
   
    %CV = im2col(YW,[bSZ bSZ],'sliding');
    %[q1,q2] = ndgrid(-20:1:20,-20:1:20);
    
    d1 = F1{t}(X(:,1),X(:,2));
    d2 = F2{t}(X(:,1),X(:,2));
    VEL = sum([d1 d2].*[d1 d2],2).^.5;
    V = [V;VEL];
    %d1 = reshape(d1,size(m1));
    %d2 = reshape(d2,size(m1));
    DELTA = oX - X;
    X = X + [d1 d2];
   
    N = sum(DELTA.*DELTA,2).^.5;
    MAG = sum(DELTA.*DELTA,2);
    N(N == 0) = 1;
    DELTA  = bsxfun(@times,DELTA,N.^-1);
    X = X + .1*bsxfun(@times,DELTA,N);
    
    toRED = VEL > VELT;
    
    
    
    xd1 = F1{t}(XM(:,1),XM(:,2));
    xd2 = F2{t}(XM(:,1),XM(:,2));
    VELM = sum([xd1 xd2].*[xd1 xd2],2).^.5;
    DELTAM = oXM - XM;
    XT = XM;
    XM = XM + [xd1 xd2];
   
    NM = sum(DELTAM.*DELTAM,2).^.5;
    MAG = sum(DELTAM.*DELTAM,2);
    
    NM(NM == 0) = 1;
    DELTAM  = bsxfun(@times,DELTAM,NM.^-1);
    XM = XM + .1*bsxfun(@times,DELTAM,NM);
    
    dT = XM - XT;
    
    PM_TMP = zeros(size(YW));
    
    XM = round(XM);
    XM((XM(:,1) <= 0),1) = 1;
    XM((XM(:,1) > size(YW,1)),1) = size(YW,1);
    XM((XM(:,2) <= 0),2) = 1;
    XM((XM(:,2) > size(YW,2)),2) = size(YW,2);
    
    
    IDX = sub2ind(size(PM_TMP),XM(:,1),XM(:,2));
    
    %diskV1 = imfilter(reshape(xd1,size(xm1)),fspecial('disk',bSZ));
    %diskV2 = imfilter(reshape(xd2,size(xm1)),fspecial('disk',bSZ));
    %diskVT = imfilter(reshape(VELM,size(xm1)),fspecial('disk',bSZ));
    %toQ = diskVT > 5;
    
    
    
    %VELM = imfilter(reshape(VELM,size(xm1)),fspecial('gaussian',[31 31],7));
    %toPROB = VELM > 5.5;
    %toPROB = toPROB(:);
    %PM_TMP = PM;
    %PM_TMP(find(toPROB)) = 1;
    %PM_TMP(IDX(toPROB)) = 1;
    %PM_TMP = imfill(PM_TMP,'holes');
    %dB = bwboundaries(logical(PM_TMP));
    
    
    
    R = X(toRED,:);
    dd = divergence(reshape(dT(:,2),size(xm1)),reshape(dT(:,1),size(xm2)));
    %dd = divergence(reshape(XM(:,2),size(xm1)),reshape(XM(:,1),size(xm2)));
    dd = imfilter(dd,fspecial('disk',bSZ),'replicate');
    %d = d + dd;
    d = dd;
    %Y = cat(3,bindVec(d),Y,Y);
    %imshow(Y,[]);
    
    tmp = ba_interp2(S(:,:,t+45),XM(:,2),XM(:,1));
    dataI(:,:,t,1) = reshape(tmp,size(xm1));
    dataI(:,:,t,2) = dd;
    dataI(:,:,t,3) = reshape(NM,size(xm1));
    %{
    tmp = permute(dataI(:,:,t,:),[1 2 4 3]);
    for k = 1:3
        tmp(:,:,k) = bindVec(tmp(:,:,k));
    end
    imshow(tmp,[]);
    %}
    
    
    %{
    figure(h1)
    imshow(d,[0 5]);
    hold on
    plot(X(:,2)-150,X(:,1)-100,'b.')
    plot(R(:,2)-150,R(:,1)-100,'r.')
    hold off
    %}
    
    %quiver(xm2(toQ),xm1(toQ),diskV2(toQ),diskV1(toQ),10)
    %   hold off
    
    %{
    figure(h2);
    O = flattenMaskOverlay(YW,logical(PM_TMP));
    imshow(O,[]);
    hold on
    %for c = 1:numel(dB)
    %    plot(dB{c}(:,2),dB{c}(:,1),'r')
    %end
    hold off
    %}
    %waitforbuttonpress
    %pause(.5)
    %hello=1
    
    t
    %drawnow
end
%%
