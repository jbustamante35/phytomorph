% read the video from disk
M = VideoReader('/mnt/tetra/nate/Drug.avi');
T = M.NumberofFrames;
H = M.Height;
W = M.Width;
S = zeros(H,W,T,'uint8');
cnt = 1;
M = VideoReader('/mnt/tetra/nate/Drug.avi');
while hasFrame(M)
    tmp = M.readFrame();
    S(:,:,cnt) = tmp(:,:,1);
    cnt = cnt + 1;
end
%% remove first N frames from stack
S(:,:,1:45) = [];
%% remove the mean frame from the stack
uS = mean(single(S),3);
sS = std(single(S),1,3);
for e = 1:size(S1,3)
    MSK = (single(S(:,:,e)) - uS) > 2*sS;
    %S1(:,:,e) = double(MSK).*double(S(:,:,e)) + double(~MSK).*uS;
    S1(:,:,e) = double(MSK).*double(S(:,:,e));
    %imshow(S1(:,:,e),[]);
    %title(num2str(e))
    %drawnow
    %imshow(S1(:,:,e),[]);
    %drawnow
    e
end
%%
%S1 = bsxfun(@minus,single(S),uS);
D1 = diff(S1,1,3);
for e = 1:size(S1,3)
    [d1,d2] = gradient(S1(:,:,e));
    dG(:,:,e) = (d1.^2 + d2.^2).^.5;
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
%%
W = [];
for e = 1:20:size(D1,3)
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = S1(:,:,e);
    H = im2colF(double(tmp),[51 51],[11 11]);
    idx = randperm(size(H,2));
    H = H(:,idx);
    H = H(:,1:5000);
    H = sort(H,1);
    H1 = H;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = D1(:,:,e);
    H = im2colF(double(tmp),[51 51],[11 11]);
    %idx = randperm(size(H,2));
    H = H(:,idx);
    H = H(:,1:5000);
    H = sort(H,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = dG(:,:,e);
    H2 = im2colF(double(tmp),[51 51],[11 11]);
    %idx = randperm(size(H,2));
    H2 = H2(:,idx);
    H2 = H2(:,1:5000);
    H2 = sort(H2,1);
    
    
    
    W = [W [H;H1;H2]];
    toc
end
%% apply beta
%% 
CHUNK = size(W,1)/3;
IDX1 = 1:CHUNK;
IDX2 = (CHUNK+1):(2*CHUNK);
IDX3 = (2*CHUNK+1):(3*CHUNK);
[U1,E1,L1] = PCA_FIT_FULL_Tws(W(IDX1,:),3);
[U2,E2,L2] = PCA_FIT_FULL_Tws(W(IDX2,:),3);
[U3,E3,L3] = PCA_FIT_FULL_Tws(W(IDX3,:),3);
L1 = L1(1);
L2 = L2(1);
L3 = L3(1);
%%
Sc1 = PCA_REPROJ_T(W(IDX1,:),E1,U1);
Sc2 = PCA_REPROJ_T(W(IDX2,:),E2,U2);
Sc3 = PCA_REPROJ_T(W(IDX3,:),E3,U3);
Sc1 = bsxfun(@times,Sc1,L1^-1);
Sc2 = bsxfun(@times,Sc2,L2^-1);
Sc3 = bsxfun(@times,Sc3,L3^-1);
%%
GMModel = fitgmdist([Sc1;Sc2]',2);
[idxL] = cluster(GMModel,[Sc1;Sc2]');
UQ = unique(idxL);
BKIDX = mode(idxL);
selIDX = idxL ~= BKIDX;
GMModel2 = fitgmdist([Sc1(:,selIDX);Sc2(:,selIDX)]',3);
%%
plot(sum(W(IDX1,:),1),Sc1(1,:),'.');
plot(Sc1(1,:),Sc1(2,:),'.');
%%
ds = 1;
cnt = 1;
vec1 = [];
for e = 1:3:size(D1,3)
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = S1(:,:,e);
    reOLD = size(tmp);
    H = im2colF(double(tmp),[51 51],[ds ds]);
    H = sort(H,1);
    Sc1 = PCA_REPROJ_T(H,E1,U1);
    Sc1 = bsxfun(@times,Sc1,L1.^-.1);
    vec1(cnt,:) = mean(Sc1,2)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = dG(:,:,e);
    reOLD = size(tmp);
    H = im2colF(double(tmp),[51 51],[ds ds]);
    H = sort(H,1);
    Sc3 = PCA_REPROJ_T(H,E3,U3);
    Sc3 = bsxfun(@times,Sc3,L3.^-.1);
    vec2(cnt,:) = mean(Sc3,2)';
    cnt = cnt + 1;
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
   
    %CV = im2col(YW,[21 21],'sliding');
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
    
    %diskV1 = imfilter(reshape(xd1,size(xm1)),fspecial('disk',21));
    %diskV2 = imfilter(reshape(xd2,size(xm1)),fspecial('disk',21));
    %diskVT = imfilter(reshape(VELM,size(xm1)),fspecial('disk',21));
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
    dd = imfilter(dd,fspecial('disk',21),'replicate');
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
