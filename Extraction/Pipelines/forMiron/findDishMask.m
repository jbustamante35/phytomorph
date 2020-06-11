function [maskImage,experimentFrame,maskStruct,background] = findDishMask(videoFile,skipN,sigSmoothValue,per,edgeCloseValue,edgeDilateValue,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skipN = 10;
    % sigSmoothValue = 11
    % per = .4;
    % edgeCloseValue = 31;
    % edgeCloseValue = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: finding pipette (rough cut).\n']);tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read every skipN frames
    [tmp,frames,FPS] = readSkipped(videoFile,-inf,skipN,inf);
    % take abs value of difference over time
    dTmp = abs(diff(tmp,1,3));
    % restack data and sum
    szT = size(dTmp);
    dTmp = reshape(dTmp,[prod(szT(1:2)) szT(3)]);
    % average change in pixel intensity for each frame
    tmpT = mean(dTmp,1);
    % smooth the abs(diff)
    tmpTs = imfilter(tmpT,fspecial('average',[1 sigSmoothValue]),'replicate');
    % sort the values
    stmpTs = sort(tmpTs);
    % get mean of the bottom/top percent
    baseLevel = mean(stmpTs(1:per*skipN));
    % get the variance of the bottom/top percent
    baseLevel_var = std(stmpTs(1:per*skipN));
    % get the normalized sigal
    nsig = (tmpTs - baseLevel)/baseLevel_var;
    pipetteSig = nsig > 10;
    pipetteSig = bwlarge(pipetteSig);
    % find the pipette frames
    pidx = find(pipetteSig);
    % get the frme to start the experiment
    experimentFrame = round(frames(pidx(end))*FPS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['end:' num2str(toc) ':finding pipette (rough cut).\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: background (rough cut).\n']);tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make min projection for background
    frames = [experimentFrame 10 Inf];
    %func = @(acc,cur).1*max(cat(3,acc,cur),[],3) + .9*min(cat(3,acc,cur),[],3);
    func = @(acc,cur)min(cat(3,acc,cur),[],3);
    minBackground = frameFunc(videoFile,frames,func,false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make max projection for background
    frames = [experimentFrame 10 Inf];
    %func = @(acc,cur).1*max(cat(3,acc,cur),[],3) + .9*min(cat(3,acc,cur),[],3);
    func = @(acc,cur)max(cat(3,acc,cur),[],3);
    maxBackground = frameFunc(videoFile,frames,func,false);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make max projection for background
    frames = [experimentFrame 10 Inf];
    N = 20;
    %func = @(acc,cur).1*max(cat(3,acc,cur),[],3) + .9*min(cat(3,acc,cur),[],3);
    func = @(acc,cur)Nmin_background(acc,cur,N);
    Nth_min = frameFunc(videoFile,frames,func,false);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make background at percent above min
    backgroundPercent = 05;
    deltaBackground = maxBackground - minBackground;
    thresholdBackground = backgroundPercent*deltaBackground + minBackground;
    bottomSignal = frameFunc(videoFile,frames,func,false);
    fprintf(['end:' num2str(toc) ':background (rough cut).\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: finding petri plate.\n']);tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the experiment slice
    [tmp,frames,FPS] = readSkipped(videoFile,experimentFrame,skipN,inf);
    % get the edges for the experiment
    ED = zeros(size(tmp));
    edgeAreaFilterValue = 300;

    %{
    % might try average then edge rather than edge then average
    for f = 1:size(tmp,3)
        [tmpJ,thJ] = edge(tmp(:,:,f),'Canny');
        tmpJ = bwareaopen(tmpJ,edgeAreaFilterValue);
        ED(:,:,f) = tmpJ;
    end
    ED = mean(double(ED),3);
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% analysis of single image rather than average of multiple analysis
    % mean of image stack - has objects near the petri plate edge
    toEdge = mean(tmp,3);
    toEdge = minBackground;
    %[d1,d2] = gradient(toEdge);
    %toEdge = (d1.^2 + d2.^2).^.5;
    [tmpJ,cannyThreshold] = edge(toEdge,'canny',[.001 .5]);
    tmpJ = bwareaopen(tmpJ,edgeAreaFilterValue);
    ED = tmpJ;


    % get the average edge mask
    uED = bindVec(ED);
    % threshold the edge mask
    uED = uED > graythresh(uED);
    % close near-by edges
    uED = imclose(uED,strel('disk',edgeCloseValue,0));
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % processing the data
    % dilate the edges
    %uED = imdilate(uED,strel('disk',edgeDilateValue,0));
    % get largest object - fail - background might be in two piece
    %uED = bwlarge(uED);
    % get only internel - fail - background might not touch edge of image
    %uED_internal = imclearborder(uED);
    %uED = uED & ~uED_internal;
    %}
    % hough transform for circle finding
    [centers,radii,metric] = imfindcircles(uED,[350 450],'EdgeThreshold',.9,'ObjectPolarity','dark','Sensitivity',.999);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['end:' num2str(toc) ':finding petri plate.\n']);
    
    if disp
        imshow(tmp(:,:,1),[]);
        %imshow(uED,[]);
        viscircles(centers(1,:), radii(1),'EdgeColor','r');
        drawnow
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ready the output
    CEN = centers(1,:);
    RAD = radii(1);
    FF = tmp(:,:,1);
     
    maskImage = zeros(size(uED));
    maskImage(round(CEN(2)),round(CEN(1))) = 1;
    maskImage = bwdist(maskImage) < RAD(1);
    
    maskStruct.center = CEN;
    maskStruct.radius = RAD;
    maskStruct.testImage = FF;
    maskStruct.maskImage = maskImage;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read upto rough cut experiment frame
    fprintf(['starting: finding pipette (fine cut).\n']);tic
    frames = [1 1 experimentFrame];
    colorMask = repmat(maskImage,[1 1 3]);
    MN = sum(maskImage(:));
    func = @(acc,cur)[acc ...
        [mean((cur(:) - minBackground(:)));...
         sum((cur(:).*maskImage(:) - minBackground(:).*maskImage(:)))/MN]];
    deltaUpto = frameFunc(videoFile,frames,func,false,[]);
    sig = bsxfun(@minus,deltaUpto,min(deltaUpto,[],2))';
    sig = mean(sig,2);
    fprintf(['end:' num2str(toc) ':finding pipette (fine cut).\n']);   
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the start of the experiment
    fprintf(['starting: finding pipette (fine cut).\n']);tic
    frames = [1 1 experimentFrame];
    func = @(acc,cur)[acc sum(sum(pipetteDetection(maskImage,minBackground,cur,6000)))>0];
    deltaUpto = frameFunc(videoFile,frames,func,false,[]);
    sig = bwlarge(deltaUpto);
    sig = imdilate(sig,strel('disk',5));
    sidx = find(sig);
    experimentFrame = sidx(end);
    fprintf(['end:' num2str(toc) ':finding pipette  (fine cut).\n']);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make min projection for background
    fprintf(['starting: full background detection.\n']);tic
    frames = [experimentFrame 10 Inf];
    func = @(acc,cur)min(cat(3,acc,imfilter(cur,fspecial('average',[3 3]))),[],3);
    %func = @(acc,cur).1*max(cat(3,acc,cur),[],3) + .9*min(cat(3,acc,cur),[],3);
    minBackground = frameFunc(videoFile,frames,func,false);
    fprintf(['end:' num2str(toc) ':full background detection.\n']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make max projection for background
    frames = [experimentFrame 10 Inf];
    %func = @(acc,cur).1*max(cat(3,acc,cur),[],3) + .9*min(cat(3,acc,cur),[],3);
    func = @(acc,cur)max(cat(3,acc,cur),[],3);
    maxBackground = frameFunc(videoFile,frames,func,false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make background at percent above min
    deltaBackground = maxBackground - minBackground;
    thresholdBackground = .05*deltaBackground + minBackground;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make max projection for background
    frames = [experimentFrame 10 Inf];
    N = 20;
    %func = @(acc,cur).1*max(cat(3,acc,cur),[],3) + .9*min(cat(3,acc,cur),[],3);
    func = @(acc,cur)Nmin_background(acc,cur,N);
    Nth_min = frameFunc(videoFile,frames,func,false);
   
    
    % old value
    background = minBackground;
    background = thresholdBackground;
    background = Nth_min(:,:,end);
    %{
    % view the pipette
    [TMP,szComplex] = readSlice(videoFile,experimentFrame,20);
    TMP = squeeze(TMP);
    for e = 1:size(TMP,3)
        imshow(TMP(:,:,e),[]);
        drawnow
        if e == 21
            waitforbuttonpress
        end
    end
    %}
    %}

    if disp
        imshow(tmp(:,:,1),[]);
        imshow(ED,[]);
        viscircles(centers(1,:), radii(1),'EdgeColor','r');
        drawnow
    end
    
end