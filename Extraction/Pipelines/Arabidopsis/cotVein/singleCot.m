function [] = singleCot(b,fileName,I,M,cropBox,sen,spurAmount,snipAmount,skeletonDilateAmount,totalMaskErodeAmount,nsz)
        
        [pth,nm,ext] = fileparts(fileName);
        
        subI = imcrop(I,cropBox);
        oI = subI;
        subI = imfilter(subI,fspecial('gaussian',[31 31],3));


    


        para.scales.value = 1;
        para.resize.value = 1;
        tmp = surKur(subI,para);
        subI = tmp(:,:,2);
        subI = imcomplement(subI);
        subI = bindVec(subI);






        %subI = imfilter(subI,fspecial('gaussian',[31 31],3));




        subM = imcrop(M,cropBox);
        subM = imfill(subM,'holes');
        
        %{
        mag = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%
        % first attempt at thresholding
        cidx = find(subM==1);
        sample = subI(cidx);
        th = graythresh(sample);
        th = th + mag*th;
        sidx = find(sample > th);
        tmpM = zeros(size(subM));
        tmpM(cidx(sidx)) = 1;
        %}

        %%%%%%%%%%%%%%%%%%%%%%%%
        % second attempt at thresholding
        T = adaptthresh(subI,sen,'ForegroundPolarity','bright','NeighborhoodSize',nsz);
        tmpM = imbinarize(subI,T);
        



       



        %%%%%%%%%%%%%%%%%%%%%%%%
        % process the cot mask
        tmpM = bwlarge(tmpM,3);
        tmpM = ~tmpM;
        tmpM = ~bwareaopen(tmpM,100);



        %tmpM = imfill(tmpM,'holes');

        veinMask = tmpM;
        [skeleton,pointSets] = extractCotVeins(veinMask,spurAmount,snipAmount);







        Cskeleton = imdilate(skeleton,strel('disk',skeletonDilateAmount));


        [skeleton,pointSets] = extractCotVeins(Cskeleton,spurAmount,snipAmount);

        %{
        cLoop = ~skeleton;
        cLoop = cLoop.*subM;
        cLoop = imerode(cLoop,strel('disk',2,0));
        %}

        [cLoopLabel,rgbC] = extractCotRegions(subM,Cskeleton,totalMaskErodeAmount);
      


        out = flattenMaskOverlay(oI,logical(subM),.3,'r');
        out = flattenMaskOverlay(out,logical(tmpM),.3,'b');
        out = flattenMaskOverlay(out,logical(skeleton),.3,'y');
        oRGB = cat(3,oI,oI,oI);
        close all
        imshow([out rgbC oRGB],[]);
        hold on
        plot(pointSets.bPoints(:,2),pointSets.bPoints(:,1),'r.')
        plot(pointSets.ePoints(:,2),pointSets.ePoints(:,1),'b.')

        %{
        for r = 1:numel(ridx)
            tmpX = sPoints(path{ridx(r),pidx(ridx(r))},2);
            tmpY = sPoints(path{ridx(r),pidx(ridx(r))},1);
            PTH = [tmpY tmpX];
            PTH = setdiff(PTH,bPoints,'rows')
            plot(PTH(:,2),PTH(:,1),'k.')
        end
        %}

        drawnow





        snipFile = [pth filesep nm '_' num2str(b) '_' num2str(sen) '.jpg'];

        saveas(gca,snipFile);

        hold off
        %waitforbuttonpress


end