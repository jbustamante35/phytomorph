%% scan for data
FilePath = '/home/nate/RILPop/';
FileList = {};
FileExt = {'csv'};
FileList = fdig(FilePath,FileList,FileExt,1);
%% load data
fidx1 = contains(FileList,'rawData');
fidx2 = contains(FileList,'RILpop');
fidx = fidx1 & fidx2;
FileList = FileList(fidx);
func = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);

clear P
basePer = .1;
disp = false;
LEN = zeros(numel(FileList),1);
G = zeros(numel(FileList),4);
err = LEN;
for e = 1:numel(FileList)
    try
        initP = zeros(1,4);


        initP(4) = 1.2;
        initP(3) = .008;

        d = csvread(FileList{e});

        [pth,nm,ext] = fileparts(FileList{e});
        didx = strfind(pth,'--');
        pth
        tnm = pth((didx(end-1)+2):(didx(end)-1));
        fidx = strfind(tnm,'_');
        l(e) = str2num(tnm((fidx(end)+1):end));

        sd = sort(d(:,2),'descend');
        initP(2) = mean(sd(1:round(numel(sd)*basePer)));
        initP(1) = mean(d(:,1));
        toFit = @(P)func(d(:,1),P);
        delta = @(P)-log(mean(normpdf((toFit(P) - d(:,2)))));


        para = fminsearch(delta,initP);

        xpos = linspace(0,1500,1500);
        vel = func(xpos,para);

        
        G(e,:) = para;
        LEN(e) = l(e);
        err(e) = delta(para);
        
        if disp
            plot(d(:,1),d(:,2),'.','MarkerSize',1);
            hold on
            plot(xpos,vel,'r');
            axis([0 1500 0 5]);
            drawnow
            hold off
        end
    catch ME
        ME
        G(e,:) = NaN*ones(1,4);
        LEN(e) = NaN;
        err(e) = NaN;
        
    end
    e
end
%% backup
GBK = G;
LENBK = LEN;
ERBK = err;
%% restore
G = GBK;
LEN = LENBK;
err = ERBK;
%% attach length
G = [G,LEN];
%% convert
conFPS = 30^-1;
conSPH = 60*60;
conPPM = 1463^-1;
G(:,2) = G(:,2)*conFPS*conSPH*conPPM;
%% remove nan
rmidx = any(isnan(G),2);
G(rmidx,:) = [];
LEN(rmidx) = [];
err(rmidx) = [];
%% outlier remove
TF = [];
for e = 1:size(G,2)
    TF(:,e) = isoutlier(G(:,e));
end
rmidx = any(TF,2);
G(rmidx,:) = [];
LEN(rmidx) = [];
err(rmidx) = [];
%% predict linear length
close all
LENP = G(:,2)*1*24;
plot(LEN+rand(size(LEN)),LENP,'.')
%% curve constraints
%%%%%%%%%%%%
% create the parameter space
np = [50 50 50 50];
masterP = [];
for e = 1:size(G,2)
    dms{e} = linspace(min(G(:,e)),max(G(:,e)),np(e));
end
[masterP(:,:,:,:,1),masterP(:,:,:,:,2),masterP(:,:,:,:,3),masterP(:,:,:,:,4)] = ndgrid(dms{:});
mSZ = size(masterP);
masterP = reshape(masterP,[prod(mSZ(1:4)) mSZ(5)]);
%%%%%%%%%%%%
% function for vel and regr
kFunc = @(X,P)P(:,2).*(1+exp(-P(:,3).*(X-P(:,1)))).^-(P(:,4).^-1);
jFunc = @(X,P)(P(:,3).*P(:,2).*exp(-P(:,3).*(X - P(:,1)))).*(P(:,4).*(exp(-P(:,3).*(X - P(:,1))) + 1).^(P(:,4).^-1 + 1)).^-1;

%%
Xmax = 3000;
Xn = 100;
Xvalues = linspace(0,Xmax,Xn);

s = 2;
P = G(s,:);
[m] = p2m(P,Xvalues);

%%
%{
%%%%%%%%%%%%
%% create symbolic funcs - base plate for the p2m function

Xmax = 3000;
Xn = 1000;
Xvalues = linspace(0,Xmax,Xn);

s = 2;
P = G(s,:);




close all
clear vel
syms xo vf k n X g f

vel(X,xo,vf,k,n) = vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);       % velocity
regr(X,xo,vf,k,n) = diff(vel,X,1);                          % regr
peak(X,xo,vf,k,n) = diff(regr,X,1);                         % where this is zero is peak
gz(X,xo,vf,k,n) = diff(peak,X,1);                           % where this is max  == diff is zero is front



f(X) = peak(X,P(1),P(2),P(3),P(4));
g(X) = regr(X,P(1),P(2),P(3),P(4));
h(X) = gz(X,P(1),P(2),P(3),P(4));


% find the peak location and value
Y = g(Xvalues);
[~,initMax] = max(Y);
peakLocation = f == 0;
peakX = vpasolve(peakLocation,X,Xvalues(initMax));
peakV = g(peakX);

% find the start points for solution
y = vpa(f(Xvalues));
[~,maxY] = max(y);
[~,minY] = min(y);
regrTipLocation = vpasolve(h==0,X,Xvalues(maxY));
regrBaseLocation = vpasolve(h==0,X,Xvalues(minY));
regrTipValue = g(regrTipLocation);
regrBaseValue = g(regrBaseLocation);
width = regrBaseLocation - regrTipLocation;

% asymetric
tipH = peakX - regrTipLocation;
baseH = regrBaseLocation - peakX;
zoneR = baseH / tipH;

% display
close all
plot(Xvalues,Y,'k');hold on
plot(peakX,peakV,'r*');
plot(regrTipLocation,regrTipValue,'g*');
plot(regrBaseLocation,regrBaseValue,'b*');
plot([regrBaseLocation regrTipLocation],[peakV peakV],'m');
title(num2str(zoneR));
%}


%%
front = h == 0;
vpasolve(front,X)
close all
figure;
plot(h(linspace(0,1000,10000)))

%%
figure;
y = g(linspace(0,1000,1000));
plot(y,'k')
hold on
plot(peakX,peakV,'r*');


%%%%%%%%%%%%
%% slow
%%%%%%%%%%%%
% compute new constrait co-ordinates
X = linspace(0,3000,1000);
u = zeros(size(masterP,1),4);
parfor e = 1:size(masterP,1)
    tic
    kField(e,:) = jFunc(X,masterP(e,:));
    u(e,:) = computeNewParameters(kField(e,:),.8);
    toc
end
u = reshape(u,[mSZ(1:numel(np)) size(u,2)]);
szU = size(u);
%%%%%%%%%%%%
%% memory
%%%%%%%%%%%%
% compute new constrait co-ordinates
X = linspace(0,3000,1000);
kField = jFunc(X,masterP);
u = computeNewParameters(kField,.8);
u = reshape(u,[mSZ(1:numel(np)) size(u,2)]);
szU = size(u);
%%%%%%%%%%%%
%% gradients
%%%%%%%%%%%%
szU = size(u);
grd = [];
for d = 1:size(u,5)
    [grd(:,:,:,:,2,d),grd(:,:,:,:,1,d),grd(:,:,:,:,3,d),grd(:,:,:,:,4,d)] = gradient(u(:,:,:,:,d),dms{2},dms{1},dms{3},dms{4});
end
%% locate the data in the cuvilinear space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toMapX = X;
dataField = jFunc(toMapX,G);
newG = computeNewParameters(dataField,.8);
close all
plot3(newG(:,1),newG(:,2),newG(:,3),'.');
newU = mean(newG);
[newS,newC,newU2,newE,newL,newERR,newLAM] = PCA_FIT_FULL(newG,3);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map mean curvi back 
close all
[oldD,con] = invMap(newU,toMapX,@(X,P)myREGR(X,P),initP);
figure;
plot3(G(:,1),G(:,2),G(:,3),'k.');
hold on
plot3(oldD(1),oldD(2),oldD(3),'r*');
%E(:,1)
cTrace = oldD;
%% make PCA spanning space and construct function to translate - work here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XalongRoot = linspace(0,3000,100);
%% pca to reduce from 4-->3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[zG,zU,zS] = zscore(G);
[latS,latC,latU,latE,latL,latERR,latLAM] = PCA_FIT_FULL(G,3);
typicalX = [diag(latLAM).^.5;1]';
typicalX = zS; % for raw trace mapping
aff = zeros(1,size(latE,2)+1);
aff(end) = 1;
affineM = [[latE,latU'];aff];
affineMz = [[latE*latLAM.^.5,latU'];aff];
velFunc = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);
regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;
bk_lat_to_org = @(p)(mtimesx(affineM,p'))';
bk_lat_to_org_simple = @(p)(mtimesx(latE,p')'+latU);
re_org_to_lat_simple = @(p)(mtimesx(latE,'T',(p-latU)'))';
re_org_to_lat = @(p)(mtimesx(pinv(affineM),p'))';
bk_nlat_to_org = @(p)(mtimesx(affineMz,p'))';
bk_org_to_nlat = @(x)(mtimesx(inv(affineMz),p'))';
z_to_org = @(p)bsxfun(@plus,bsxfun(@times,p,zS),zU);
org_to_z = @(p)bsxfun(@times,bsxfun(@minus,p,zU),zS.^-1);
% map the data points to the curvilinear space from the parameter space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;[mG] = p2m(G,XalongRoot);toc
tic;[mGn] = p2mn(G);toc
% NOTE: which conversion to use !!!!!
[zG,curU,curSIG] = zscore(mG);
[curS,curC,curU,curE,curL,curERR,curLAM] = PCA_FIT_FULL(mG,3);
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~,curE,~,~,curLAM] = PCA_FIT_FULL(zscore(mG),3);
for e = 1:size(curE,2)
    tmp = curE(:,e)'.*curSIG;
    curE(:,e) = (tmp / norm(tmp))';
end
tmpC = PCA_REPROJ(mG,curE,curU);
curLAM = diag(std(tmpC,1,1).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
    tmp = curE(:,1);
    tmp(2) = 0;
    tmp = tmp / norm(tmp);
    curE(:,1) = tmp;
%}

%norCur = @(x)(x - curU).*curSIG.^-1;
%inorCur = @(x)(x.*curSIG)+curU;
norCur = @(x)bsxfun(@times,bsxfun(@minus,x,curU),curSIG.^-1); % one is added for affine mapping
%inorCur = @(x)(x.*curSIG(1:size(curC,2)));
curU = mean(mG);
%plot(curC(:,1),C(:,2),'.');
for e = 1:size(mG,2)
    figure;
    plot(mG(:,e),mGn(:,e),'.')
    waitforbuttonpress
    close all
end
% plot pca in curvi

%% try to inverse map from curvi to pca above - not yet working - WORKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initP = [0;0;0;1];
traceC = {};
deltaVec = [];
for pcN = 1:3
    
    
    N = 21;
    alpha = 1;
    mag = 2.7;
    %mag = 1;
    CL = {'r','g'};
    
    
    alongEigenVec = alpha*linspace(-mag*(curLAM(pcN,pcN).^.5),mag*(curLAM(pcN,pcN).^.5),N);
    %alongEigenVec = alpha*linspace(-2^.5,2^.5,N);
    pointsAlong = bsxfun(@plus,(curE(:,pcN)*alongEigenVec)',curU);
    %pointsAlong = (curE(:,pcN)*alongEigenVec)';
    pointsAlong = norCur(pointsAlong);
    
    %pointsAlong = pointsAlong(1,:);
    
    clear testP
    [testP] = m2p(pointsAlong,@(v)norm(v),initP,XalongRoot,bk_nlat_to_org,norCur);
    
    % measure error in point finding
    delta = norCur(p2mn(bk_nlat_to_org(testP))) - pointsAlong;
    deltaVec(:,pcN) = sum(delta.^2,2);
    traceC{pcN} = bk_nlat_to_org(testP);
    pcN
end

%% re-para the trace curve
for e = 1:numel(traceC)
    AtraceC{e} = arcLength(traceC{e},'spec',size(traceC{e},1));
end
%%
close all
%% plot the search methods
% note need to bkProject into the G space

close all
CL = {'r','g','b'};
for e = 1:numel(AtraceC)
    plot3(G(:,1),G(:,2),G(:,3),'k.');hold on;
    plot3(AtraceC{e}(:,1),AtraceC{e}(:,2),AtraceC{e}(:,3),CL{e},'Linewidth',2);
end
%% N choose K
close all
P = nchoosek(1:4,2);
LABEL = {'x0','vf','k','n','len'};
for e = 1:size(P,1)
    figure;
    tX = G(:,P(e,:));
    plot(tX(:,1),tX(:,2),'k.');
    hold on;
    for c = 1:numel(AtraceC)
        tCurve = AtraceC{c}(:,P(e,:));
        plot(tCurve(:,1),tCurve(:,2),CL{c},'Linewidth',2);
    end
    xlabel(LABEL{P(e,1)});
    ylabel(LABEL{P(e,2)});
    waitforbuttonpress
    close all
end
%% scan the eigen vectors
close all
X = linspace(0,2500,1000);
func = @(X,P)P(2).*(1+exp(-P(3).*(X-P(1)))).^-(P(4).^-1);
for e = 1:numel(AtraceC)
    figure
    for n = 1:3:size(AtraceC{e},1)
        Y = func(X,AtraceC{e}(n,:));
        gY = gradient(Y);
        
        conFPS = 30^-1;
        conSPH = 60*60;
        conPPM = 1463^-1;
        % convert vel
        Y = Y * conFPS * conSPH *conPPM;
        % convert regr
        gY = gY * conFPS * conSPH * 100;
        
        %G(:,2) = G(:,2)*conFPS*conSPH*conPPM;
        yyaxis left
        plot(X,Y,[CL{e} '-']);
        
        yyaxis right
        plot(X,gY,['k' '-']);
        hold on
    end
    waitforbuttonpress
end
%%

%% try to inverse map through the sub space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear P
velFunc_p = @(X,P)velFunc(X,bk_lat_to_org(P));
regrFunc_p = @(X,P)regrFunc(X,bk_lat_to_org(P));
regrPer = .8;
% this line was converted so that number of points along the domain do not
% influence the results
%toMetricSpace_p = @(P)extractREGRmetrics(@(P)regrFunc_p(XalongRoot,P),P,regrPer);
toMetricSpace_p = @(P)p2m(bk_lat_to_org(P),XalongRoot);
toMetricSpace_p2 = @(P)p2mn(bk_lat_to_org(P));
toMetricSpace_pn = @(P)p2mn(bk_nlat_to_org(P));
toMetricSpace_p3 = @(P)p2mn(bk_lat_to_org_simple(P));
toMetricSpace_raw = @(P)p2mn(P);
toMetricSpace_ZT = @(P)p2mn(z_to_org(P));
%%
close all
pTest = [0;0;0;0]';
if all(bk_lat_to_org(pTest) == [mean(G,1),1]);fprintf(['Mean from parameter space maps to mean of data space.\n']);end
%% map from metrix space to parameter space - here latent parameter space is a stand-in for the 
% TEST inverse mappping 
pTest = [0;0;0;0]';
[testP] = m2p(pTest,@(v)(norm(v)),initP,XalongRoot,bk_nlat_to_org,norCur);
norCur(p2mn(bk_nlat_to_org(testP)))
norm(pTest-norCur(p2mn(bk_nlat_to_org(testP))))
%%
regrMeanPCA = regrFunc_p(XalongRoot,pTest);
regrMeanOriginal = regrFunc(XalongRoot,mean(G,1));
plot(XalongRoot,regrMeanOriginal,'k');hold on
plot(XalongRoot,regrMeanPCA+.00001,'r');
hold on;
u = toMetricSpace_p(pTest);
figure;
plot(XalongRoot,regrMeanPCA,'k');hold on
plot(u(1),u(2),'g*');
%% TRACE METHOD
% find the co-ordinate [0,0,0,0]@metrix/curvilinear in parameter/lat space
curP = curU;
[latP] = m2p(curP,@(v)(norm(v)),initP,XalongRoot,bk_lat_to_org,norCur);
bk_lat_to_org(latP)
%%
nTrace = 50;
%nTrace = 1;
dv = curE(:,1);
alpha = .1;
%traceInLat1_pos = [0;0;0;1]';
%traceInLat1_pos = [0;0;0]';
%traceInLat1_pos = latP / latP(end);
%traceInLat1_pos = latP;
%traceInLat1_pos = re_org_to_lat(AtraceC{1}(11,:));
traceInLat1_pos = org_to_z(AtraceC{1}(11,1:size(G,2)));
[traceInLat1_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_pos,dv,nTrace,alpha,typicalX);
traceInLat1_neg = org_to_z(AtraceC{1}(11,1:size(G,2)));
[traceInLat1_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat1_neg,-dv,nTrace,alpha,typicalX);
%
nTrace = 50;
dv = curE(:,2);
alpha = .1;
traceInLat2_pos = org_to_z(AtraceC{2}(11,1:size(G,2)));
[traceInLat2_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat2_pos,dv,nTrace,alpha,typicalX);
traceInLat2_neg = org_to_z(AtraceC{2}(11,1:size(G,2)));
[traceInLat2_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat2_neg,-dv,nTrace,alpha,typicalX);
%
nTrace = 50;
dv = curE(:,3);
alpha = .1;
traceInLat3_pos = org_to_z(AtraceC{3}(11,1:size(G,2)));
[traceInLat3_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat3_pos,dv,nTrace,alpha,typicalX);
traceInLat3_neg = org_to_z(AtraceC{3}(11,1:size(G,2)));
[traceInLat3_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat3_neg,-dv,nTrace,alpha,typicalX);
%{
LABEL = {'x0','vf','k','n','len'};
CON = {'peakPosition','peak','Width','asym'};
%}
%% trace out components of grid
nTrace = 50;
DV = eye(4);
%dv = [0;0;0;1];
dv = DV(:,3);
alpha = .1;
traceInLat4_pos = org_to_z(AtraceC{3}(11,1:size(G,2)));
[traceInLat4_pos] = traceVec_ver2(toMetricSpace_ZT,traceInLat4_pos,dv,nTrace,alpha,typicalX);
traceInLat4_neg = org_to_z(AtraceC{3}(11,1:size(G,2)));
[traceInLat4_neg] = traceVec_ver2(toMetricSpace_ZT,traceInLat4_neg,-dv,nTrace,alpha,typicalX);
%%

%%
close all
plot3(G(:,1),G(:,2),G(:,3),'k.');
hold on

proj = z_to_org(traceInLat1_pos);
plot3(proj(:,1),proj(:,2),proj(:,3),'r','LineWidth',2);
proj = z_to_org(traceInLat1_neg);
plot3(proj(:,1),proj(:,2),proj(:,3),'r','LineWidth',2);
plot3(AtraceC{1}(:,1),AtraceC{1}(:,2),AtraceC{1}(:,3),'m')


proj = z_to_org(traceInLat2_pos);
plot3(proj(:,1),proj(:,2),proj(:,3),'g','LineWidth',2);
proj = z_to_org(traceInLat2_neg);
plot3(proj(:,1),proj(:,2),proj(:,3),'g','LineWidth',2);
plot3(AtraceC{2}(:,1),AtraceC{2}(:,2),AtraceC{2}(:,3),'y')

proj = z_to_org(traceInLat3_pos);
plot3(proj(:,1),proj(:,2),proj(:,3),'b','LineWidth',2);
proj = z_to_org(traceInLat3_neg);
plot3(proj(:,1),proj(:,2),proj(:,3),'b','LineWidth',2);
plot3(AtraceC{3}(:,1),AtraceC{3}(:,2),AtraceC{3}(:,3),'c')

proj = z_to_org(traceInLat4_pos);
plot3(proj(:,1),proj(:,2),proj(:,3),'k','LineWidth',2);
proj = z_to_org(traceInLat4_neg);
plot3(proj(:,1),proj(:,2),proj(:,3),'k','LineWidth',2);
plot3(AtraceC{3}(:,1),AtraceC{3}(:,2),AtraceC{3}(:,3),'c')


hold off
view([90 0])
drawnow
%%
%{
%% template for trace routine above
close all
mag = .1; % works
mag = .1; % works
%mag = 1.5;

for e = 1:200
   
    %[f,j] = extractMetricsAtP(toMetricSpace_p,traceInLat(end,:),1*ones(size(pTest)),typicalX);
    
    tic
    [f,j] = extractMetricsAtP(toMetricSpace_p2,traceInLat(end,:),1*ones(size(pTest)),typicalX);
    toc
    
    %[U,S,V] = svd(j');
    
    %U = U(1:3,1:3);
    %V = V(1:3,1:3);
    %S = S(1:3,1:3);
    %J = inv(U*S*V');
    
    jT = j';
    %jT(:,end) = [];
    % focus on the mapping - bttom is the source top is the target
    % [d{target}/d{source}] = Tij = i,target j,source
    dc = [(pinv(jT)*dv)'];
    tt = jT*dc';
    %dc = inv(j')*dv;
%    dc = J*dv
    %dc = dc / dc(end);
    norm(dc)
    dc = mag*(dc / norm(dc));
    %dc = [dc 0];
    
    %dc(1) = 0;
    %dc(end) = 0;
    %dc(2:3) = mag*(dc(2:3) / norm(dc(2:3)));
     %dc = [dc 0];
    %{
    dc = mag*(dc / norm(dc));
    dc = [dc 0];
    %}
    
    
    
    traceInLat(e+1,:)  = traceInLat(e,:) + dc;
    e
    
    
   
    %proj = PCA_BKPROJ(traceInLat(:,1:(end-1)),latE,latU);
    %{
    plot(G(:,2),G(:,3),'k.');
    hold on
    plot(proj(:,2),proj(:,3),'r')
    %}
    
    %plot(proj(:,2))
    
    %drawnow
    %hold off
    
    plot3(G(:,1),G(:,2),G(:,3),'k.');
    hold on
    proj = PCA_BKPROJ(traceInLat(:,1:(end-1)),latE,latU);
    %proj = PCA_BKPROJ(traceInLat(:,1:(end-1)),latE,0);
    plot3(proj(:,1),proj(:,2),proj(:,3),'r')
    hold off
    view([90 0])
    drawnow
    
    %verfi = p2mn(proj);
end
%}
%%

[f,j] = extractMetricsAtP(func,metricFunc,p,h);

%% make nD interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
masterP = reshape(masterP,mSZ);
for i = 1:size(grd,5)
    for j = 1:size(grd,6)
        interP{i,j} = griddedInterpolant(masterP(:,:,:,:,1),...
            masterP(:,:,:,:,2),masterP(:,:,:,:,3),masterP(:,:,:,:,4),...
            grd(:,:,:,:,i,j),'cubic');
    end
end
%%

close all
alpha = 1;
beta = .2;
alongN = 1000;

bvec = newE(:,1);
n = 51;
dl = linspace(0,newLAM(1,1).^.5,alongN);
[cTracePos1] = traceVec(interP,grd,dms,cTrace,bvec,dl,alpha,beta);
dl = linspace(0,-newLAM(1,1).^.5,alongN);
[cTraceNeg1] = traceVec(interP,grd,dms,cTrace,bvec,dl,alpha,beta);
cTracePosT1 = [flip(cTraceNeg1,1);cTracePos1];
scTracePosT1 = imfilter(cTracePosT1,ones(n,1)/n,'replicate');


bvec = newE(:,2);
dl = linspace(0,newLAM(2,2).^.5,alongN);
[cTracePos2] = traceVec(interP,grd,dms,cTrace,bvec,dl,alpha,beta);
dl = linspace(0,-newLAM(2,2).^.5,alongN);
[cTraceNeg2] = traceVec(interP,grd,dms,cTrace,bvec,dl,alpha,beta);
cTracePosT2 = [flip(cTraceNeg2,1);cTracePos2];
scTracePosT2 = imfilter(cTracePosT2,ones(n,1)/n,'replicate');

%{
for l = 1:numel(dl)
    for m = 1:size(grd,5)
        for p = 1:size(grd,6)
            dM(m,p) = interpn(dms{1},dms{2},dms{3},dms{4},grd(:,:,:,:,m,p),...
                pStack(l,1),pStack(l,2),pStack(l,3),pStack(l,4),'cubic');
        end
    end
    vec = (inv(dM)'*newE(:,1)*dl(l)*alpha)';
    n = norm(vec);
    if n~= 0
        vec = vec / norm(vec);
    end
    pStack(l+1,:) = pStack(l,:) + beta*vec;
end
%}
close all
figure
plot3(G(:,1),G(:,2),G(:,3),'k.');
hold on
plot3(oldD(1),oldD(2),oldD(3),'r*');
drawnow

plot3(scTracePosT1(:,1),scTracePosT1(:,2),scTracePosT1(:,3),'r')
plot3(scTracePosT2(:,1),scTracePosT2(:,2),scTracePosT2(:,3),'g')
%%
plot3(cTracePos1(:,1),cTracePos1(:,2),cTracePos1(:,3),'r')
plot3(cTraceNeg1(:,1),cTraceNeg1(:,2),cTraceNeg1(:,3),'r')

plot3(cTracePos2(:,1),cTracePos2(:,2),cTracePos2(:,3),'g')
plot3(cTraceNeg2(:,1),cTraceNeg2(:,2),cTraceNeg2(:,3),'g')
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
X = linspace(0,2000,2000);
kField = jFunc(X,masterP);
u = computeNewParameters(kField,.8);
u = reshape(u,[mSZ(1:4) size(u,2)]);
szU = size(u);
for d = 1:size(u,5)
    d
    [grd(:,:,:,:,1,d),grd(:,:,:,:,2,d),grd(:,:,:,:,3,d),grd(:,:,:,:,4,d)] = gradient(u(:,:,:,:,d));
end

%% locate the data in the cuvilinear space
toMapX = X;%linspace(0,3000,5000);
dataField = jFunc(toMapX,G);
newG = computeNewParameters(dataField,.8);
close all
plot3(newG(:,1),newG(:,2),newG(:,3),'.');
[newS newC newU newE newL newERR newLAM] = PCA_FIT_FULL(newG,3);

%% test procedure for mapping mean to mean
%{ 
figure;
plot3(newG(:,1),newG(:,2),newG(:,3),'.');
initP = mean(G,1);

testF = jFunc(toMapX,initP);
testP = computeNewParameters(testF,.8);
xp = invMap(testP,X,@(X,P)myREGR(X,P),initP);


testF1 = jFunc(toMapX,xp);
testP2 = computeNewParameters(testF1,.8);
xp2 = invMap(testP2,X,@(X,P)myREGR(X,P),initP);
close all

figure;
plot(myREGR(X,initP)+.001,'k.')
hold on
plot(myREGR(X,xp)-.001,'r')
plot(myREGR(X,xp2),'b')
waitforbuttonpress
close all
%}

%% map mean curvi back 
close all
[oldD,con] = invMap(newU,toMapX,@(X,P)myREGR(X,P),initP);
figure;
plot3(G(:,1),G(:,2),G(:,3),'k.');
hold on
plot3(oldD(1),oldD(2),oldD(3),'r*');
%%
curvi1 = newE(:,2)*linspace(-newLAM(1,1).^.5,newLAM(1,1).^.5,30);
curvi1 = bsxfun(@plus,curvi1,newU')';
hold on
plot3(curvi1(:,1),curvi1(:,2),curvi1(:,3),'b.');hold on;
%waitforbuttonpress
[zG,uG,sigG] = zscore(G);
[UU,con] = invMap(curvi1,toMapX,@(X,P)myREGR(X,P),initP);
%%
%{
r = figure;
plot(myREGR(X,UU))
hold on
plot(myREGR(X,mean(G,1)),'r')
plot(myREGR(X,initP),'g.')
waitforbuttonpress
close(r)
%}

%{
u = reshape(u,[prod(szU(1:4)) szU(end)]);
for e = 1:size(curvi1,1)
    tmpD = bsxfun(@minus,u,curvi1(e,:));
    tmpD = sum(tmpD.^2,2);
    [~,midx(e)] = min(tmpD);
end
u = reshape(u,szU);
clear curvi1P
[curvi1P(:,1),curvi1P(:,2),curvi1P(:,3),curvi1P(:,4)] = ind2sub(szU(1:(end-1)),midx);
%}


%dataFieldTest = jFunc(toMapX,curvi1P);
%newGTest = computeNewParameters(dataFieldTest,.8);
curvi1P = UU;
figure;
plot3(G(:,1),G(:,2),G(:,3),'k.')
hold on
%plot3(UU(1),UU(2),UU(3),'r*')
%plot3(dms{1}(curvi1P(:,1)),dms{2}(curvi1P(:,2)),dms{3}(curvi1P(:,3)),'r')
plot3((curvi1P(:,1)),(curvi1P(:,2)),(curvi1P(:,3)),'r')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slices of u
close all
P = nchoosek(1:size(G,2),2);
LABEL = {'x0','vf','k','n','len'};
CON = {'peakPosition','peak','Width','asym'};
szU = size(u);
singleMode = true;
multiMode = 1;
DIMS = [3 4];

for e = 1:size(P,1)
    pv = [P(e,:) setdiff(1:(ndims(u)-1),P(e,:))];
    pv = [pv ndims(u)];
    ud = permute(u,pv);
    upSample = 2;
    
    
    if ~singleMode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % multipage
        for l = 1:size(u,DIMS(multiMode))

            % single snap slice
            plot(G(:,P(e,2)),G(:,P(e,1)),'k.')
            xlabel(LABEL{P(e,2)})
            ylabel(LABEL{P(e,1)})

            title(CON{m});
            hold on
            CL = {'r','g','b','k'};
            for m = 1:szU(end)
                
                
                if multiMode == 1
                    slice = ud(:,:,l,25,m);
                else
                    slice = ud(:,:,25,l,m);
                end
                
                
                sliceX = dms{P(e,1)};
                sliceY = dms{P(e,2)};

                sliceXi = interp1(linspace(0,1,numel(sliceX)),sliceX,linspace(0,1,upSample*numel(sliceX)));
                sliceYi = interp1(linspace(0,1,numel(sliceY)),sliceY,linspace(0,1,upSample*numel(sliceY)));
                [sliceY,sliceX] = ndgrid(sliceX,sliceY);
                [sliceYi,sliceXi] = ndgrid(sliceXi,sliceYi);
                slice = interp2(sliceX,sliceY,slice,sliceXi,sliceYi);
                contour(sliceXi,sliceYi,slice,10,CL{m});
                
            end
            if multiMode == 1
                title({['R->PeakLocation:G->PeakValue:B->Peak Width'],...
                    [LABEL{pv(3)} '<-C:' LABEL{pv(4)} '<-Sweep']});
            else
                title({['R->PeakLocation:G->PeakValue:B->Peak Width'],...
                    [LABEL{pv(4)} '<-C:' LABEL{pv(3)} '<-Sweep']});
            end
            
            drawnow
            hold off
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    if singleMode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % single snap slice
        plot(G(:,P(e,2)),G(:,P(e,1)),'k.')
        xlabel(LABEL{P(e,2)})
        ylabel(LABEL{P(e,1)})

        hold on
        CL = {'r','g','b','k'};
        for m = 1:szU(end)

            %title(CON{m});
            slice = ud(:,:,25,25,m);
            sliceX = dms{P(e,1)};
            sliceY = dms{P(e,2)};

            sliceXi = interp1(linspace(0,1,numel(sliceX)),sliceX,linspace(0,1,upSample*numel(sliceX)));
            sliceYi = interp1(linspace(0,1,numel(sliceY)),sliceY,linspace(0,1,upSample*numel(sliceY)));
            [sliceY,sliceX] = ndgrid(sliceX,sliceY);
            [sliceYi,sliceXi] = ndgrid(sliceXi,sliceYi);
            slice = interp2(sliceX,sliceY,slice,sliceXi,sliceYi);
            contour(sliceXi,sliceYi,slice,10,CL{m});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
        hold off
        title(['R->PeakLocation:G->PeakValue:B->Peak Width']);
    end
    waitforbuttonpress
    %u = ipermute(u,pv);
    %slice = u(:,:,25,25);
end
%% scatter plot PC vs data
close all
[Z,uean,sig] = zscore(G);
[S C U E L ERR LAM]  = PCA_FIT_FULL(Z,3);

LABEL = {'x0','vf','k','n','len'};

P = nchoosek(1:size(G,2),2);
close all
figure;
LX = linspace(-100,100,100);
for e = 1:size(P,1)
   
    figure;
    vx = LABEL{1,P(e,1)};
    vy = LABEL{1,P(e,2)};
    tX = G(:,P(e,1));
    tY = G(:,P(e,2));
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make prob map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpX = tX;
    tmpY = tY;
    NX = 1000;
    NY = 1000;
    rangeX = linspace(min(tX(:)),max(tX(:)),NX);
    rangeY = linspace(min(tY(:)),max(tY(:)),NY);
    N = [NX NY];
    numL = 8;
    samN = 100;
    for loop = 1:100
        
        
        [pM] = makeProbMap(tmpX,tmpY,N,50);
        
        imshow(pM,[]);
        drawnow
        
        newX = [];
        newY = [];
        for n = 1:samN
            [newX(n),newY(n)] = pinky(rangeX,rangeY,pM);
        end
        hold on
        plot(bindVec(newX(:))*size(pM,2),bindVec(newY(:))*size(pM,1),'r.')
       
        minZ = min(pM(:)); 
        maxZ = max(pM(:));
        zStep = linspace(minZ,maxZ,numL);
        C = contourc(rangeX,rangeY,pM,zStep);
        [dP] = extractLevelsFromContourMatrix(C);
        maxArea = 0;
        
        for c = 1:numel(dP)
            cX = dP{c}(1,:)';
            cY = dP{c}(2,:)';
            cX = [cX;cX(1)];
            cY = [cY;cY(1)];
            
            cX = nbindVec(cX,min(tX),max(tX))*size(pM,2);
            cY = nbindVec(cY,min(tY),max(tY))*size(pM,1);
            
            BW = poly2mask(cX,cY,size(pM,2),size(pM,1));
            
            area(c) = sum(BW(:));
            
            if area(c) > maxArea
                maxArea = area(c);
                maxContour = [cX,cY];
            end
            plot(cX,cY,'r');
            hold on
            
            
        end
        
        
        % analyze skeleton
        BW = poly2mask(maxContour(:,1),maxContour(:,2),size(pM,2),size(pM,1));
        SKEL = bwmorph(BW,'skeleton',inf);
        [SK(:,2),SK(:,1)] = find(SKEL);
        branchPoints = bwmorph(SKEL,'branchPoints');
        endPoints = bwmorph(SKEL,'endPoints');
        [BR(:,2),BR(:,2)] = find(branchPoints);
        [EP(:,2),EP(:,1)] = find(endPoints);
        dT = pdist2(SK,SK);
        
        
        
        plot(maxContour(:,1),maxContour(:,2),'b');
        plot(sk1,sk2,'y.');
        
        
        
        
        drawnow
        hold off
        %waitforbuttonpress
        tmpX = [tX;newX'];
        tmpY = [tY;newY'];
        loop
    end
    
    
    
    hold on
    plot(samx,samy,'.')
    
    
    
    probMapDomain = probMapDist < 35;
    probMapDomain = imfill(probMapDomain,'holes');
    probMapDomain = bwlarge(probMapDomain);
    
    
    % smooth the domain
    probMapDomain = imfilter(probMapDomain,fspecial('gaussian',31,7),'replicate');
    probMapDomain = probMapDomain > .6;
    probMapDomain = imfill(probMapDomain,'holes');
    probMapDomain = logical(probMapDomain);
    probMapDist = (-bwdist(probMapDomain) + bwdist(~probMapDomain));
    for t = 1:200
        probMapDist = probMapDist - .01*del2(probMapDist);
    end
    %figure;
    %imshow(probMapDist,[]);
    %waitforbuttonpress
    probMapDomain = probMapDist > 20;
    probMapDomain = bwlarge(probMapDomain);
    probMapDomainSkel = bwmorph(probMapDomain,'skeleton',inf);
    probMapDomainSkel = bwmorph(probMapDomainSkel,'spur',300);
    [sy,sx] = find(probMapDomainSkel==1);
    sx = bindVec(sx);
    sy = bindVec(sy);
    sx = sx * (max(tX(:)) - min(tX(:))) + min(tX(:));
    sy = sy * (max(tY(:)) - min(tY(:))) + min(tY(:));
    %imshow(probMapDomain,[]);
    %waitforbuttonpress
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    hold on
    plot(tX,tY,'.')
    CL = {'r','g','b'};
    %tU = uean([P(e,1) P(e,2)]);
    %plot(tU(2),tU(1),'m*')
    for v = 1:3
        vec = E(P(e,:),v)';
        vec = vec / norm(vec);
        vec = bsxfun(@times,LX'*vec,sig(P(e,:)));
        vec = bsxfun(@plus,vec,uean(P(e,:)));
        plot(vec(:,1),vec(:,2),CL{v});
        hold on
    end
    %plot(sx,sy,'m.')
    xlabel(vx);
    ylabel(vy);
    %contour(probMap,10)
    axis([min(tX) max(tX) min(tY) max(tY)]);
    waitforbuttonpress
end
close all
%% predicted length of root
close all
yLEN = LEN - mean(LEN);
W = yLEN\C;
preLEN = C*W';
plot(yLEN+mean(LEN),preLEN+mean(LEN),'k.');
hold on
plot([0 6],[0 6],'r')
%%  scatter scores
for e = 1:size(C,2)
    figure;
    Y = LEN;
    X = C(:,e);
    [xcor,pval] = corr(X,Y);
    Y = Y + rand(size(Y));
    plot(C(:,e),Y,'.');
    ylabel('init L')
    xlabel(['pc' num2str(e)]);
    title(['cor-' num2str(xcor) ':pval-' num2str(pval)])
    waitforbuttonpress
    close all
end

for e = 1:size(G,2)
    figure;
    Y = LEN;
    Y = Y + rand(size(Y));
    X = G(:,e);
    [xcor,pval] = corr(X,Y);
    plot(X,Y,'.');
    ylabel('init L')
    xlabel(LABEL{e});
    title(['cor-' num2str(xcor) ':pval-' num2str(pval)])
    waitforbuttonpress
    close all
end
%% subplots for sweeping eigna vectors
func = @(X,xo,vf,k,n)vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);
close all
[sweepD] = sweepPCA(C,E,U,diag(LAM),1:3,30);


XY = linspace(0,2000,20000);

h1 = figure;
%h2 = figure;
close all

for c = 1:size(sweepD,1)
    
    %h1 = figure;
    %h2 = figure;
    tmp = squeeze(sweepD(c,:,:));
    
    tmp = bsxfun(@times,tmp,sig);
    tmp = bsxfun(@plus,tmp,uean);
    p = [];
    v = [];
    
    
    
    for n = 1:size(tmp,1)
        %Y = func(XY,wholeD(s,1),tmp(n,1),tmp(n,2),tmp(n,3));
        Y = func(XY,tmp(n,1),tmp(n,2),tmp(n,3),tmp(n,4));
        velP(c,n) = Y(end);
        %Y = gradient(Y);
        %figure(h1);
        %plot(XY,Y);
        %hold all
        
        Y = gradient(Y);
        
        [v(c,n),p(c,n)] = max(Y);
        
        msk = Y > v(c,n)*.80;
        
        
        loc = find(msk);
        
        H1(c,n) = p(c,n) - loc(1);
        H2(c,n) = loc(end) - p(c,n);
        
        
        
        
        WID(c,n) = max(XY(loc)) - min(XY(loc));
        
        
        p(c,n) = XY(p(c,n));
        %figure(h2)
        %plot(XY,Y);
        %hold all
        
    end
    
    
    [MM,alpha] = max(v(c,:));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for eigen vector one - plot cell lengths
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if c == 1
        para = tmp(alpha,:);
        g = find(v(c,:) > (min(v(c,:)) + .8*(max(v(c,:)) - min(v(c,:)))));
        q1 = tmp(g(1),:);
        q2 = tmp(g(end),:);

        initLength = .005;
        initPos = .005;
        mag = 1;
        
        initPos = 100*1463^-1;
        initLength = .020;
        
        initP = [mag*initPos mag*(initPos)+initLength];
        
        dt = 1;
        TAU = 24*3;

        
        nT = round(TAU / dt);
        tmD = linspace(0,TAU,nT)/24;
        tmD = tmD';
        
        initP = [initPos initLength];
        
        [Lcenter] = genCellLength(para,initP,dt,TAU);
        [Lless] = genCellLength(q1,initP,dt,TAU);
        [Lgreater] = genCellLength(q2,initP,dt,TAU);

        %Lcenter = Lcenter(1:numel(tmD));
        %Lless = Lless(1:numel(tmD));
        %Lgreater = Lgreater(1:numel(tmD));
        
        %Lcenter = Lcenter / diff(initP);
        %Lless = Lless / diff(initP);
        %Lgreater = Lgreater / diff(initP);

        h1 = figure;
        plot(Lless(:,1),Lless(:,2),'r');hold on
        plot(Lcenter(:,1),Lcenter(:,2),'g');
        plot(Lgreater(:,1),Lgreater(:,2),'b');
        waitforbuttonpress
        close(h1);
    end
    
    
    
    x = p(c,:)/1463;
    x = (x - min(x)) / max(x);
    
    reg = v(c,:);
    reg = (reg - min(reg)) / max(reg);
    
    
    baseP = (c - 1)*3;
    
    
    subplot(3,3,baseP + 1)
    %figure;
    plot(x*100,reg*100);
    hold on
    if c == 1
        plot(x(g(1))*100,reg(g(1))*100,'r*');
        plot(x(alpha)*100,reg(alpha)*100,'g*');
        plot(x(g(end))*100,reg(g(end))*100,'b*');
    end
    
    title('REGR peak location vs REGR peak value');
    if c == 1
        xlabel('Percent change in peak position');
        ylabel('Percent change in peak value');
    end
    
    
    subplot(3,3,baseP + 2)
    y = velP(c,:);
    y = (y - min(y)) / max(y);
    plot(x*100,y*100)
    title('REGR peak location vs VEL value');
    if c == 1
        xlabel('Percent change in peak position');
        ylabel('Percent change in velocity');
    end
    
    
    subplot(3,3,baseP + 3)
    w = WID(c,:);
    h1 = H1(c,:);
    h2 = H2(c,:);
    w = (w - min(w)) / max(w);
    h1 = (h1 - min(h1)) / max(h1);
    h2 = (h2 - min(h2)) / max(h2);
    CL = {'r','g','b'};
    
    plot(x*100,w*100,CL{1});
    hold on
    plot(x*100,h1*100,CL{2});
    plot(x*100,h2*100,CL{3});
    
    
    title('REGR peak location vs WIDTH');
    if c == 1
        xlabel('Percent change in peak position');
        ylabel('Percent change in width');
    end
    
    
    
    
    LEG{1} = '80 percent width';
    LEG{2} = 'tip half';
    LEG{3} = 'base half';
    %legend(LEG);
    
    
    waitforbuttonpress
    
    %waitforbuttonpress
    %close all
    %figure(h2);
    
    
end
%%
D = readtext('~/Cvi_parameters.csv');
wholeD = cell2mat(D(2:end,3:6));
subD = cell2mat(D(2:end,4:6));
Z = zscore(subD);
[Z,uean,sig] = zscore(wholeD);
%%
close all
func = @(X,xo,vf,k,n)vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);
exY = func(XY,wholeD(1,1),wholeD(1,2),wholeD(1,3),wholeD(1,4));
plot(XY,exY)
evaY = func(wholeD(1,1),wholeD(1,1),wholeD(1,2),wholeD(1,3),wholeD(1,4));
b1 = evaY^wholeD(1,4)
b2 = (wholeD(1,2)^wholeD(1,4))/2
%%

toP = nchoosek(1:4,2);
toP=[flip(toP,2);toP];
LAB = {'Xo','Vf','k','n'};
XY = linspace(0,2000,100);
func = @(X,xo,vf,k,n)vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);
oPath = '/home/nate/forSpalding/';
%toP = [2 3];

toA = wholeD;
toA = Z;

X = toA;
close all
plot3(X(:,1),X(:,2),X(:,3),'k.');
totC = 3;
[S C U E L ERR LAM] = PCA_FIT_FULL(X,totC);


for e = 1:size(toP,1)
    
    
    close all
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    %h2 = figure;
    
    tX = toA(:,toP(e,1));
    tY = toA(:,toP(e,2));
    tP = toA;
    
    b = robustfit(tX,tY);
    lX = linspace(.5*min(tX),1.5*max(tX),100);
    [stX,sidx] = sort(tX);
    stY = tY(sidx);
    
    
    tP = tP(sidx,:);
    
    tCur = [];
    for p = 1:size(tP,1)
        tCur(p,:)= func(XY,tP(p,1),tP(p,2),tP(p,3),tP(p,4));
    end
    
    
    
    figure(h1);
    subplot(1,2,1)
    hold on
    SKIP = 10;
    IDX = 1:SKIP:size(tCur,1);
    LEG = {};
    CL = {'r','g','b','k'};
    for p = 1:numel(IDX)
        plot(XY,tCur(IDX(p),:),CL{p})
        LEG{p} = num2str(stX(IDX(p)));
    end
    
    
    legend(LEG)
    %waitforbuttonpress
    figure(h1);
    subplot(1,2,2)
    preY = lX*b(2) + b(1);
    plot(tX,tY,'.');
    hold all
    xlabel(LAB{toP(e,1)});
    ylabel(LAB{toP(e,2)});
    plot(lX,preY,'r');
    
    
    
    
    CL = {'m','g','c'};
    hold on
    
    
    for v = 1:size(E,2)
        
        
        mXX = linspace(-1000,1000,3);
        be = E(toP(e,:),v);
        %be = be.*sig(toP(e,:))';
        be = be/norm(be');
        myG = (be*mXX)';
        %myG = myG + (uean(toP(e,:)));
        plot(myG(:,1),myG(:,2),CL{v});
        
        
        % for kicks
        if v == 1 && e == 4
            pause(2)
            figure;
            bv = E(:,v);
            pc1 = toA*bv;
            sweepV = std(pc1);
            L1 = linspace(-sweepV,sweepV,5);
            %vec = bsxfun(@plus,bv*L1,uean');
            vec = bv*L1;
            regr = [];
            p = [];
            for i = 1:size(vec,2)
                p(:,i) = vec(:,i);
                p(:,i) = p(:,i).*sig' + uean';
                regr(i,:) = func(XY,p(1,i),p(2,i),p(3,i),p(4,i));
                %regr(i,:) = gradient(regr(i,:));
                plot(XY,regr(i,:),'r')
                hold on
            end
            
            bv = E(:,v);
            pc1 = toA*bv;
            
            toRe = 1:4;
            toRe = setdiff(toRe,toP(e,:));
            bv(toRe) = 0;
            vec = bv*L1;
            regr = [];
            for i = 1:size(vec,2)
                pq = vec(:,i);
                pq = pq.*sig' + uean';
                
                regr(i,:) = func(XY,pq(1),pq(2),pq(3),pq(4));
                %regr(i,:) = gradient(regr(i,:));
                plot(XY,regr(i,:),'b')
                hold on
            end
            
            
        end
        
    end
    
    
    
    
    
    %plot(uean(toP(e,1)),uean(toP(e,2)),'m*')
    plot(0,0,'m*')
    axis([lX(1),lX(end),min(tY),max(tY)])
    CL = {'ro','go','bo','ko'};
    for p = 1:numel(IDX)
        plot(stX(IDX(p)),stY(IDX(p)),CL{p})
    end
    
    
    
    
    hold off
    waitforbuttonpress
    fileName = [oPath LAB{toP(e,1)} '-' LAB{toP(e,2)} '.tif'];
    saveas(gca,fileName);
    
end
%% sweep x0 hold n 
close all

XY = linspace(0,2000,100);

nsubD = [];
func = @(X,xo,vf,k,n)vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);
ux0 = mean(wholeD(:,1));
ux0 = .6*ux0;


LB = .8*min(wholeD(:,2:4),[],1);
UB = 1.2*max(wholeD(:,2:4),[],1);
PSS = 2000;
for e = 1:size(wholeD,1)
    
    
    
    
    UtmpY = func(XY,ux0,wholeD(e,2),wholeD(e,3),wholeD(e,4));
    tmpY = func(XY,wholeD(e,1),wholeD(e,2),wholeD(e,3),wholeD(e,4));
    
    
    %pltF =  @(P)plotT(tmpY,func,P);
    %ux0 = wholeD(e,1);
    toMin = @(P)norm(tmpY - func(XY,ux0,P(1),P(2),P(3)));
    %toMin = @(P)max(tmpY - func(XY,ux0,P(1),P(2),P(3)));
    
    
    options = optimoptions(@particleswarm,'Display','iter','SwarmSize',PSS);
    
    %x = fminsearch(toMin,wholeD(e,2:4),options);
    
    %x = fminunc(toMin,wholeD(e,2:4));%,options);
    
    x = particleswarm(toMin,3,LB,UB,options);
    
    
    options = optimset('Display','iter','TolX',10^-1000,'TolFun',10^-1000);
    
    x = fminsearch(toMin,x,options);
    
    newY = func(XY,ux0,x(1),x(2),x(3));
    
    plot(XY,tmpY,'k')
    hold on
    plot(XY,UtmpY,'r')
    plot(XY,newY,'b')
    
    OLD = wholeD(e,2:4);
    NEW = x;
    
    %title([num2str(OLD) '--' num2str(NEW)])
    title(num2str(e))
    hold off
    drawnow
    %waitforbuttonpress
    nsubD(e,:) = x;
end

%%
X = subD;
X = nsubD;
X = wholeD;
X = Z;
close all
plot3(X(:,1),X(:,2),X(:,3),'k.');
totC = 3;
[S C U E L ERR LAM] = PCA_FIT_FULL(X,totC);
hold on
CL = {'r' , 'g' ,'b'};
scale = diag(LAM);
%%
for loop = 1:4
    ROT = linspace(0,360,100);
    cnt = 1;

    for i = 1:size(X,1)

        plot3(X(:,1),X(:,2),X(:,3),'k.');
        plot3(U(1),U(2),U(3),'k*')
        hold on
        for e = 1:totC
            e
            tmp = [U;U+2*scale(e)*E(:,e)'];
            plot3(tmp(:,1),tmp(:,2),tmp(:,3),CL{e})
        end
        plot3(U(1),U(2),U(3),'c')

        P = (E*C')' + U;

        plot3(P(i,1),P(i,2),P(i,3),'b.')

        dis = P - X;

        %for e = 1:size(P,1)
            tmp = [P(i,:);X(i,:)];
            plot3(tmp(:,1),tmp(:,2),tmp(:,3),'g')

            tmp = [U;P(i,:)];
            plot3(tmp(:,1),tmp(:,2),tmp(:,3),'b')
        %end
        %waitforbuttonpress
       % hold off
        axis
        cnt = 1;
        for r = 1:numel(ROT)
            view(32,ROT(cnt));
            cnt = cnt + 1;
            drawnow;
            pause(.2);
        end
        
    end
end
%% look at SIM vs real
close all
for s = 1:size(wholeD,1)
   % Ysim = func(XY,wholeD(s,1),S(s,1),S(s,2),S(s,3));
    tmpPara = S(s,:);
    tmpPara = (sig.*tmpPara + uean);
    Ysim = func(XY,tmpPara(1,1),tmpPara(1,2),tmpPara(1,3),tmpPara(1,4));
    Yreal = func(XY,wholeD(s,1),wholeD(s,2),wholeD(s,3),wholeD(s,4));
    plot(XY,Yreal,'k');
    hold on
    plot(XY,Ysim,'r');
    axis([0 max(XY) 0 4])
    hold off
    waitforbuttonpress
end
%%

close all
[sweepD] = sweepPCA(C,E,U,diag(LAM).^.5,1:3,30);


XY = linspace(0,2000,20000);

h1 = figure;
%h2 = figure;
close all

for c = 1:size(sweepD,1)
    
    %h1 = figure;
    %h2 = figure;
    tmp = squeeze(sweepD(c,:,:));
    
    tmp = bsxfun(@times,tmp,sig);
    tmp = bsxfun(@plus,tmp,uean);
    p = [];
    v = [];
    
    
    
    for n = 1:size(tmp,1)
        %Y = func(XY,wholeD(s,1),tmp(n,1),tmp(n,2),tmp(n,3));
        Y = func(XY,tmp(n,1),tmp(n,2),tmp(n,3),tmp(n,4));
        velP(c,n) = Y(end);
        %Y = gradient(Y);
        %figure(h1);
        %plot(XY,Y);
        %hold all
        
        Y = gradient(Y);
        
        [v(c,n),p(c,n)] = max(Y);
        
        msk = Y > v(c,n)*.80;
        
        
        loc = find(msk);
        
        H1(c,n) = p(c,n) - loc(1);
        H2(c,n) = loc(end) - p(c,n);
        
        
        
        
        WID(c,n) = max(XY(loc)) - min(XY(loc));
        
        
        p(c,n) = XY(p(c,n));
        %figure(h2)
        %plot(XY,Y);
        %hold all
        
    end
    
    
    [MM,alpha] = max(v(c,:));
    
    para = tmp(alpha,:);
    g = find(v(c,:) > (min(v(c,:)) + .8*(max(v(c,:)) - min(v(c,:)))));
    q1 = tmp(g(1),:);
    q2 = tmp(g(end),:);
    
    
    initP = [100 110];
    dt = 1;
    TAU = 2000;
    
    [Lcenter] = genCellLength(para,initP,dt,TAU);
    [Lless] = genCellLength(q1,initP,dt,TAU);
    [Lgreater] = genCellLength(q2,initP,dt,TAU);
    
    h1 = figure;
    plot(Lless,'r');hold on
    plot(Lcenter,'g*');
    plot(Lgreater,'b');
    waitforbuttonpress
    close(h1);
    
   
    
    
    x = p(c,:)/1463;
    x = (x - min(x)) / max(x);
    
    reg = v(c,:);
    reg = (reg - min(reg)) / max(reg);
    
    
    baseP = (c - 1)*3;
    subplot(3,3,baseP + 1)
    %figure;
    plot(x*100,reg*100);
    hold on
    plot(x(g(1))*100,reg(g(1))*100,'r*');
    plot(x(alpha)*100,reg(alpha)*100,'g*');
    plot(x(g(end))*100,reg(g(end))*100,'b*');
    
    title('REGR peak location vs REGR peak value');
    xlabel('Percent change in peak position');
    
    
    subplot(3,3,baseP + 2)
    y = velP(c,:);
    y = (y - min(y)) / max(y);
    plot(x*100,y*100)
    title('REGR peak location vs VEL value');
    xlabel('Percent change in peak position');
    ylabel('Percent change in velocity');
    
    
    subplot(3,3,baseP + 3)
    w = WID(c,:);
    h1 = H1(c,:);
    h2 = H2(c,:);
    w = (w - min(w)) / max(w);
    h1 = (h1 - min(h1)) / max(h1);
    h2 = (h2 - min(h2)) / max(h2);
    CL = {'r','g','b'};
    
    plot(x*100,w*100,CL{1});
    hold on
    plot(x*100,h1*100,CL{2});
    plot(x*100,h2*100,CL{3});
    
    
    title('REGR peak location vs WIDTH');
    xlabel('Percent change in peak position');
    ylabel('Percent change in width');
    LEG{1} = '80 percent width';
    LEG{2} = 'tip half';
    LEG{3} = 'base half';
    %legend(LEG);
    
    
    waitforbuttonpress
    
    %waitforbuttonpress
    %close all
    %figure(h2);
    
    
end
%%
P = nchoosek(1:4,2);

close all
figure;
LX = linspace(-100,100,100);
for e = 1:size(P,1)
   
    figure;
    vx = D{1,2+P(e,1)};
    vy = D{1,2+P(e,2)};
    tX = wholeD(:,P(e,1));
    tY = wholeD(:,P(e,2));
    hold on
    plot(tX,tY,'.')
    CL = {'r','g','b'};
    tU = uean([P(e,1) P(e,2)]);
    plot(tU(2),tU(1),'m*')
    for v = 1:3
        vec = E(P(e,:),v)';
        vec = vec / norm(vec);
        vec = bsxfun(@times,LX'*vec,sig(P(e,:)));
        vec = bsxfun(@plus,vec,uean(P(e,:)));
        plot(vec(:,1),vec(:,2),CL{v});
        hold on
    end
    xlabel(vx);
    ylabel(vy);
    axis([min(tX) max(tX) min(tY) max(tY)]);
    %waitforbuttonpress
end
%%

P = nchoosek(1:3,2);
P = [flip(P,2);P];

for e = 1:size(P,1)
    vx = D{1,3+P(e,1)};
    vy = D{1,3+P(e,2)};
    tX = subD(:,P(e,1));
    tY = subD(:,P(e,2));
    b = robustfit(tX,tY);
    l = linspace(min(tX),max(tX));
    y = b(2)*l+b(1);
    plot(tX,tY,'k.');
    hold on
    plot(l,y,'r');
    xlabel(vx)
    ylabel(vy)
    waitforbuttonpress
   close all
end
%%
func=@(X,xo,vf,k,n)vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);
%%
XY = linspace(0,2000,100);

%%
s = 1;
Y = func(XY,wholeD(s,1),wholeD(s,2),wholeD(s,3),wholeD(s,4));
tmpS = S(s,:);
tmpS = tmpS.*sig + uean;
Ysim = func(XY,tmpS(s,1),tmpS(s,2),tmpS(s,3),tmpS(s,4));
close all
plot(XY,Y,'k');
hold all
plot(XY,Ysim,'r.');