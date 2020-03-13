fileName = '/mnt/snapper/nate/ART_test/Exp149-Substack_49-621.tif';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and stack data
% can be either a CZI of TIF file
filt = fspecial('gaussian',[41 41],11);
STACK = czi_tif_loader(fileName,filt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STACK = double(STACK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize the raw stack
[SIG] = normalizeRawStack(STACK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make whole plant mask
SS = std(STACK,1,3);
% make the mean signal of the image
U = mean(STACK,3);
% make the plant mask based on the standard dev
plantMask = bindVec(SS) > graythresh(bindVec(SS))*.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the gradient data and filter
filt_SIZE = [11 11];
filt_STD = 3;
[d1,d2,d3] = gradient(SIG);
for tm = 1:size(SIG,3)
    d1(:,:,tm) = imfilter(d1(:,:,tm),fspecial('gaussian',filt_SIZE,filt_STD),'replicate');
    d2(:,:,tm) = imfilter(d2(:,:,tm),fspecial('gaussian',filt_SIZE,filt_STD),'replicate');
end
dN = (d1.^2 + d2.^2).^-.5;
d1 = d1.*dN;
d2 = d2.*dN;
clear d3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(SIG);
SIG = reshape(SIG,[prod(sz(1:2)) sz(3)]);
d1 = reshape(d1,[prod(sz(1:2)) sz(3)]);
d2 = reshape(d2,[prod(sz(1:2)) sz(3)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oSIG = reshape(SIG,sz);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DN = 5;
[s1,s2] = ndgrid(1:size(plantMask,1),1:size(plantMask,2));
ds = mod(s1,DN) == 0 & mod(s2,DN) == 0;
didx = find(ds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plantMask = bwlarge(plantMask);
pidx = find(plantMask);
pidx = intersect(pidx,didx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = .3;
disp = false;
[n2,n1] = ndgrid(linspace(-20,20,50),linspace(-20,20,50));
domainSZ = size(n1);
domain = [n1(:) n2(:) ones(size(n1(:))) ones(size(n2(:)))]';
TMI = 150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q1 = zeros([domainSZ TMI numel(pidx)]);
%Q2 = zeros([domainSZ TMI numel(pidx)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subSIG = SIG(pidx,:);
subd1 = d1(pidx,:);
subd2 = d2(pidx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subP = [];
[subP(:,1),subP(:,2)] = ind2sub(sz(1:2),pidx);
tmg = false;
%% watch oSIG
for e = 1:size(oSIG,3)
    imshow(oSIG(:,:,e),[0 1]);
    title(num2str(e))
   
    drawnow
end
%%%%%%%%%%%%
%% view a Transformation - JUNK VIEW
[w1,w2,V] = impixel(STACK(:,:,30));

threshold = .2;
e = sub2ind(size(plantMask),w2,w1);
sigStack = [SIG(e,:);d1(e,:);d2(e,:)];
[riseDomain,recoverDomain] = generateBinaryDomain(sigStack,threshold);
p = [w2 w1];
%riseDomain = imdilate(riseDomain,strel('disk',31,0));
[rise_affineSequence] = generateAffineSequence(sigStack,p,riseDomain);
tm = find(riseDomain);

close all
for t = 1:numel(tm)
    g = rise_affineSequence(1:2,1,t);
    n = rise_affineSequence(1:2,2,t);
    m = rise_affineSequence(1:2,4,t);
    
    
    sim = zeros(size(plantMask));
    sim(w2-20,w1-20) = 1;
    sim = double(bwdist(sim));
    
    
    [l1,l2] = gradient(oSIG(:,:,tm(t)));
    %[l1,l2] = gradient(sim);
    %{
    imshow(l1,[-.2 .2]);
    waitforbuttonpress
    %}
    l1 = ba_interp2(l1,w1,w2);
    l2 = ba_interp2(l2,w1,w2);
    NC = norm([l1 l2])
    l1 = l1 / NC;
    l2 = l2 / NC;
    
    test(t) = ba_interp2(oSIG(:,:,tm(t)),w1,w2);
    l1
    imshow(oSIG(:,:,tm(t)),[0 1]);
    %imshow(sim,[ ]);
    hold on
    quiver(w1,w2,g(1),g(2),20);
    %quiver(w1,w2,n(1),n(2),20);
    quiver(w1,w2,l1,l2,20,'g');
    plot(w1,w2,'go');
    plot(m(1),m(2),'r.');
    title(num2str(tm(t)))
    drawnow
    hold off
end
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

kp = [];
parfor e = 1:numel(pidx)
    try
        TTtm = clock;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        sigStack = [subSIG(e,:);subd1(e,:);subd2(e,:)];
        etm = etime(clock,tm);
        if tmg;fprintf(['stack data signals :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        %[p(1),p(2)] = ind2sub(sz(1:2),pidx(e));
        p = subP(e,:);
        etm = etime(clock,tm);
        if tmg;fprintf(['translate p-idx to p-xy :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        [riseDomain,recoverDomain] = generateBinaryDomain(sigStack,threshold);
        etm = etime(clock,tm);
        if tmg;fprintf(['generate binary domains  :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        [rise_affineSequence] = generateAffineSequence(sigStack,p,riseDomain);
        %[recover_affineSequence] = generateAffineSequence(sigStack,p,recoverDomain);
        etm = etime(clock,tm);
        if tmg;fprintf(['generate affine sequences :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_affineSequence] = rescaleData(rise_affineSequence,TMI);
        %[recover_affineSequence] = rescaleData(recover_affineSequence,TMI);
        etm = etime(clock,tm);
        if tmg;fprintf(['rescale affine data :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        riseFrames = oSIG(:,:,(riseDomain));
        recoverFrames = oSIG(:,:,(recoverDomain));
        etm = etime(clock,tm);
        fprintf(['index frames2 :' num2str(etm) '\n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_data_out] = coreSample(oSIG,domain,rise_affineSequence,domainSZ,disp);
        %[recover_data_out] = coreSample(oSIG,domain,recover_affineSequence,domainSZ,disp);
        etm = etime(clock,tm);
        if tmg;fprintf(['sample core data :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_data_out] = rescaleData(rise_data_out,200);
        [recover_data_out] = rescaleData(recover_data_out,200);
        etm = etime(clock,tm);
        fprintf(['rescale core data :' num2str(etm) '\n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}

        Q1(:,:,:,e) = rise_data_out;
        %Q2(:,:,:,e) = recover_data_out;
        etm = etime(clock,TTtm);
        fprintf(['total time :' num2str(etm) '\n'])
        kp(e) = true;
    catch
        kp(e) = false;
    end
end
%% threshold wave
close all
nQ1 = zeros(size(tmpQ1));
for e = 1:size(tmpQ1,4)
    for t = 1:size(Q1,3)
        msk = tmpQ1(:,:,t,e) > tmpQ1(25,25,t,e);
        %msk = bwdist(msk) - bwdist(~msk);
        %imshow(msk,[]);
        %drawnow
        nQ1(:,:,t,e) = msk;
    end
    e
    size(tmpQ1,4)
end
%% make and remove level 1
tmpQ1 = Q1(:,:,:,find(kp));
qSZ = size(tmpQ1);
reduce1 = 3;
tmpQ1 = reshape(tmpQ1,[prod(qSZ(1:2)) prod(qSZ(3:4))]);
[U_rise1,E_rise1,L_rise1] = PCA_FIT_FULL_Tws(tmpQ1,reduce1);
C_rise1 = PCA_REPROJ_T(tmpQ1,E_rise1,U_rise1);
C_rise1 = reshape(C_rise1,[reduce1 qSZ(3:4)]);
qSZ1 = size(C_rise1);
reduce2 = 3;
tmpC1 = reshape(C_rise1,[prod(qSZ1(1:2)) prod(qSZ1(3))]);
[U_rise2,E_rise2,L_rise2] = PCA_FIT_FULL_Tws(tmpC1,reduce2);
C_rise2 = PCA_REPROJ_T(tmpC1,E_rise2,U_rise2);
C_rise2 = reshape(C_rise2,[reduce2 qSZ1(3)]);
TF = isoutlier(C_rise2(1,:));
tmpQ1 = reshape(tmpQ1,qSZ);
tmpQ1(:,:,:,TF) = [];

%% rinse and repeat
close all
qSZ = size(tmpQ1);
reduce1 = 3;
tmpQ1 = reshape(tmpQ1,[prod(qSZ(1:2)) prod(qSZ(3:4))]);
[U_rise1,E_rise1,L_rise1] = PCA_FIT_FULL_Tws(tmpQ1,reduce1);
C_rise1 = PCA_REPROJ_T(tmpQ1,E_rise1,U_rise1);
C_rise1 = reshape(C_rise1,[reduce1 qSZ(3:4)]);
qSZ1 = size(C_rise1);
reduce2 = 3;
tmpC1 = reshape(C_rise1,[prod(qSZ1(1:2)) prod(qSZ1(3))]);
[U_rise2,E_rise2,L_rise2] = PCA_FIT_FULL_Tws(tmpC1,reduce2);
C_rise2 = PCA_REPROJ_T(tmpC1,E_rise2,U_rise2);
C_rise2 = reshape(C_rise2,[reduce2 qSZ1(3)]);
TF = isoutlier(C_rise2(1,:));
tmpQ1 = reshape(tmpQ1,qSZ);
tmpQ1(:,:,:,TF) = [];
plot3(C_rise2(1,:),C_rise2(2,:),C_rise2(3,:),'.')
%% mask only - try one pass
tmpQ1 = nQ1;
close all
qSZ = size(tmpQ1);
reduce1 = 3;
tmpQ1 = reshape(tmpQ1,[prod(qSZ(1:2)) prod(qSZ(3:4))]);
[U_rise1,E_rise1,L_rise1] = PCA_FIT_FULL_Tws(tmpQ1,reduce1);
C_rise1 = PCA_REPROJ_T(tmpQ1,E_rise1,U_rise1);
C_rise1 = reshape(C_rise1,[reduce1 qSZ(3:4)]);
qSZ1 = size(C_rise1);
reduce2 = 3;
tmpC1 = reshape(C_rise1,[prod(qSZ1(1:2)) prod(qSZ1(3))]);
[U_rise2,E_rise2,L_rise2] = PCA_FIT_FULL_Tws(tmpC1,reduce2);
C_rise2 = PCA_REPROJ_T(tmpC1,E_rise2,U_rise2);
C_rise2 = reshape(C_rise2,[reduce2 qSZ1(3)]);
TF = isoutlier(C_rise2(1,:));
tmpQ1 = reshape(tmpQ1,qSZ);
tmpQ1(:,:,:,TF) = [];
plot3(C_rise2(1,:),C_rise2(2,:),C_rise2(3,:),'.')
%%
close all
for p = 1:3
    [sweepD] = sweepPCA([0 0 0],E_rise1,U_rise1',3*L_rise1.^.5,p,10);
    sweepD = squeeze(reshape(sweepD,[size(sweepD,1) size(sweepD,2) 50 50]));
    for e = 1:size(sweepD,1)
        imshow(squeeze(sweepD(e,:,:)),[]);
        drawnow
        waitforbuttonpress
    end
end
%%
close all
CL = {'r' 'g' 'b'};
for p = 1:3
    plot(C_rise1(p,:,100)',CL{p})
    hold on
end
%%
close all
LEG = {'a1' 'a2' 'a3'};
for p = 1:3
    [sweepD] = sweepPCA([0 0 0],E_rise2,U_rise2',3*L_rise2.^.5,p,10);
    sweepD = squeeze(reshape(sweepD,[size(sweepD,1) size(sweepD,2) 3 150]));
    for e = 1:size(sweepD,1)
        for k = 1:3
            plot(squeeze(sweepD(e,k,:))',CL{k});
            hold on
        end
        if e == 1
            legend(LEG)
        end
           legend(LEG)
        title(num2str(p))
        hold on
        waitforbuttonpress
    end
    hold off
end
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plantMask = bwlarge(plantMask);
%[~,BOX] = imcrop(plantMask);
%[~,BOX] = imcrop(oSIG(:,:,100));
%plantMask_C = zeros(size(plantMask));
%plantMask_C(BOX(2):(BOX(2)+BOX(4)),BOX(1):BOX(1)+BOX(3)) = 1;
pidx = find(plantMask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subSIG = SIG(pidx,:);
subd1 = d1(pidx,:);
subd2 = d2(pidx,:);
subP = [];
[subP(:,1),subP(:,2)] = ind2sub(sz(1:2),pidx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = [];
tmg = false;
kp = [];
close all
parfor e = 1:numel(pidx)
    try
        TTtm = clock;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        sigStack = [subSIG(e,:);subd1(e,:);subd2(e,:)];
        etm = etime(clock,tm);
        if tmg;fprintf(['stack data signals :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        %[p(1),p(2)] = ind2sub(sz(1:2),pidx(e));
        p = subP(e,:);
        etm = etime(clock,tm);
        if tmg;fprintf(['translate p-idx to p-xy :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        [riseDomain,recoverDomain] = generateBinaryDomain(sigStack,threshold);
        etm = etime(clock,tm);
        if tmg;fprintf(['generate binary domains  :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        [rise_affineSequence] = generateAffineSequence(sigStack,p,riseDomain);
        %[recover_affineSequence] = generateAffineSequence(sigStack,p,recoverDomain);
        etm = etime(clock,tm);
        if tmg;fprintf(['generate affine sequences :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_affineSequence] = rescaleData(rise_affineSequence,TMI);
        %[recover_affineSequence] = rescaleData(recover_affineSequence,TMI);
        etm = etime(clock,tm);
        if tmg;fprintf(['rescale affine data :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm  = clock;
        riseFrames = oSIG(:,:,(riseDomain));
        recoverFrames = oSIG(:,:,(recoverDomain));
        etm = etime(clock,tm);
        fprintf(['index frames2 :' num2str(etm) '\n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_data_out] = coreSample(oSIG,domain,rise_affineSequence,domainSZ,disp);
        %[recover_data_out] = coreSample(oSIG,domain,recover_affineSequence,domainSZ,disp);
        etm = etime(clock,tm);
        if tmg;fprintf(['sample core data :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_data_out] = ThresholdcoreSample(rise_data_out);
        %[recover_data_out] = coreSample(oSIG,domain,recover_affineSequence,domainSZ,disp);
        etm = etime(clock,tm);
        if tmg;fprintf(['threshold sample core data :' num2str(etm) '\n']);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        [rise_data_out] = rescaleData(rise_data_out,200);
        [recover_data_out] = rescaleData(recover_data_out,200);
        etm = etime(clock,tm);
        fprintf(['rescale core data :' num2str(etm) '\n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}

        
        
        C(:,e) = funky(rise_data_out,U_rise1,E_rise1,U_rise2,E_rise2);

        
        etm = etime(clock,TTtm);
        fprintf(['total time :' num2str(etm*numel(pidx)/60/60/12) '\n'])
        kp(e) = true;
    catch ME
        getReport(ME)
        kp(e) = false;
    end
end
%% make image
close all
imageFun = zeros([size(plantMask) size(C,1)]);
sz = size(imageFun);
imageFun = reshape(imageFun,[prod(sz(1:2)) sz(3)]);
imageFun(pidx,:) = C';
imageFun = reshape(imageFun,sz);

d_plantMask = imerode(plantMask,strel('disk',21,0));
d_plantMask = plantMask;
for k = 1:size(imageFun,3)
    tmp = imageFun(:,:,k);
    %TF = isoutlier(tmp);
    %tmp(TF) = 0;
    nimageFun(:,:,k) = bindVec(d_plantMask.*tmp);
end


imshow(nimageFun,[]);
for k = 1:size(imageFun,3)
    figure;
    imshow(nimageFun(:,:,k),[]);
end
%%
SKIP = 10;
plot3(C(1,1:SKIP:end),C(2,1:SKIP:end),C(3,1:SKIP:end),'.')
GMModel = fitgmdist(C',3);
%%
close all

kidx = GMModel.cluster(C');
LAB = zeros([size(plantMask)]);
LAB(pidx) = kidx;
%%
LAB = label2rgb(LAB);
imshow(LAB,[]);
%%
CL = {'r' 'g' 'b'};
close all
kidx = GMModel.cluster(C');
LAB = zeros([size(plantMask)]);
LAB(pidx) = kidx;
for t = 1:size(oSIG,3)
    a = oSIG(:,:,t);
    tmpM = a < 0;
    a(tmpM) = 0;
    %imshow(tmpM,[]);
    %drawnow
    
    out = plantMask.*bindVec(plantMask.*a);
    for k = 1:3
        out = flattenMaskOverlay(out,LAB==k,.1,CL{k});
    end
    imshow(out,[]);
    
end


