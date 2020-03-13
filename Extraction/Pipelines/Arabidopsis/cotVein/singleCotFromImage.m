function [dataPoint] = singleCotFromImage(dataPoint,sen,gapClose,lengthFilter,smallHoleFilter,snipAmount,skeletonDilateAmount,totalMaskErodeAmount,nsz)
    try
        %[pth,nm,ext] = fileparts(fileName);
        
        %subI = imcrop(I,cropBox);
        % make the use image
        subI = dataPoint.alignedImage.*dataPoint.alignedMask;
        subM = dataPoint.alignedMask;
        
        oI = subI;
        subI = imfilter(subI,fspecial('gaussian',[31 31],3));


        subM = imfill(subM,'holes');
        eD = 15;
        blendEdgeD = imerode(subM,strel('disk',eD,0));
        blendEdge = imfilter(blendEdgeD,fspecial('gaussian',[31 31],8));
        
        
        surKurSmoothValue = 1;
        para.scales.value = surKurSmoothValue;
        para.resize.value = 1;
        tmp = surKur(subI,para);
        
        %f1 = imfilter(tmp(:,:,1),fspecial('gaussian',[31 31],3));
        %f2 = imfilter(tmp(:,:,2),fspecial('gaussian',[31 31],3));
        
       
        subI = tmp(:,:,2);
        subI = imcomplement(subI);
        subI = bindVec(subI);
        %subI = imfilter(subI,fspecial('gaussian',[31 31],3));
     
        %(bindVec(blendEdge.*subI)),[])
        %subI = bindVec(blendEdge.*subI);


        %subM = imcrop(M,cropBox);
       
        
        %{
        mag = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % first attempt at thresholding
        cidx = find(subM==1);
        sample = subI(cidx);
        th = graythresh(sample);
        th = th + mag*th;
        sidx = find(sample > th);
        tmpM = zeros(size(subM));
        tmpM(cidx(sidx)) = 1;
        %}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % second attempt at thresholding
        T = adaptthresh(subI,sen,'ForegroundPolarity','bright','NeighborhoodSize',nsz);
        tmpM = imbinarize(subI,T);
        

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%
        % try to remove things that are too close the edge
        % problem is that the bottom of the cut is edge too
        %%%%%%%%%%%%%%%%%%%%%%%%
        ridx = [];
        bwD = bwdist(~subM);
        Re = regionprops(logical(tmpM),'PixelIdxList');
        for e = 1:numel(Re)
            vec = bwD(Re(e).PixelIdxList);
            if min(vec) < 5
                ridx = [ridx e];
            end
        end
        Re(ridx) = [];
        tmpM = zeros(size(tmpM));
        for e = 1:numel(Re)
            tmpM(Re(e).PixelIdxList) = 1;
        end
        tmpM = logical(tmpM);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % knock out edge objects by erosion of edge toward center
        tmpM = blendEdgeD.*tmpM;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % process the cot mask
        % keep the largest two
        tmpM = bwlarge(tmpM,2);
        tmpM = bwareaopen(tmpM,600);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init the first vein mask and extract the skeleton
        veinMask = tmpM;
        [skeleton,pointSets] = extractCotVeins(veinMask,gapClose,lengthFilter,smallHoleFilter,snipAmount);
        % only use the largest object
        skeleton = bwlarge(skeleton,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % skeleton extraction
        doFlag = true;
        while doFlag
            % store the skeleton for comparing 
            preS = skeleton;
            % extract the skeleton as new-skeleton
            [newSkeleton,pointSets,Adj] = extractCotVeins(skeleton,gapClose,lengthFilter,smallHoleFilter,snipAmount);
            % post skeleton is the one post processing
            postS = newSkeleton;
            
            % compare the pre to the post
            [r,c] = find(abs(preS - postS) ~= 0 );

            %{
            DI = cat(3,preS,postS,abs(preS - postS));
            imshow(DI,[]);
            hold on
            plot(c,r,'r.')
            drawnow
            waitforbuttonpress
            %}
            if all(newSkeleton(:) == skeleton(:))
                doFlag = false;
                skeleton = newSkeleton;
            else
                skeleton = newSkeleton;
            end
        end

        %{
        cLoop = ~skeleton;
        cLoop = cLoop.*subM;
        cLoop = imerode(cLoop,strel('disk',2,0));
        %}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dilate the skeleton for extraction of the holes etc
        Cskeleton = imdilate(skeleton,strel('disk',skeletonDilateAmount));
        % extract the cot regions
        [cLoopLabel,rgbC,holesImage] = extractCotRegions(subM,Cskeleton,totalMaskErodeAmount);
      




        inHole = [];
        for e = 1:size(pointSets.ePoints,1)
            inHole(e) = holesImage(pointSets.ePoints(e,1),pointSets.ePoints(e,2));
        end

        if any(inHole == 1)
            hidx = find(inHole);
            %{
            for h = 1:numel(hidx)
                plot(pointSets.ePoints(hidx(h),2),pointSets.ePoints(hidx(h),1),'r*')
            end
            %}
            for h = 1:numel(hidx)
                [path,pathcost] = searchAlongNetwork(Adj,pointSets.sPoints,pointSets.bPoints,pointSets.ePoints(hidx(h),:));
                [~,toRemove] = min(pathcost);
                rPath = path{toRemove};
                skeleton = removeBranch(skeleton,pointSets.sPoints,rPath,pointSets.bPoints);
            end
            [skeleton,pointSets,Adj] = extractCotVeins(skeleton,gapClose,lengthFilter,smallHoleFilter,snipAmount);
  
            Cskeleton = imdilate(skeleton,strel('disk',skeletonDilateAmount));
            [cLoopLabel,rgbC,holesImage]= extractCotRegions(subM,Cskeleton,totalMaskErodeAmount);
        end

        
            
       
        [numberHoles] = countHoles(holesImage);


        dataPoint.labelImage = rgbC;
        dataPoint.holesImage = holesImage;
        dataPoint.numberHoles = numberHoles;


        dataPoint.veinNetwork = ...
            cotGraph(Adj,pointSets.sPoints,pointSets.ePoints,pointSets.bPoints,skeleton);

        %snipFile = [pth filesep nm '_' num2str(b) '_' num2str(sen) '.jpg'];

        %saveas(gca,snipFile);

        %hold off
        %waitforbuttonpress
    catch ME
        ME
    end

end