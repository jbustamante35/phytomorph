function [cLoopLabel,rgbC,holesImage] = extractCotRegions(subM,skeleton,totalMaskErodeAmount)
    try
        cLoop = ~skeleton;
        useM = bwlarge(subM);
        for bor = 1:4
            useM(:,1) = 0;
            useM = imrotate(useM,90);
        end
        useM = imerode(useM,strel('disk',totalMaskErodeAmount,0));
        cLoop = cLoop.*useM;

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % make the holes image
        % find the first pixel and remove the object 
        % associated with it
        % find is in raster order
        holesImage = ~cLoop;
        fidx = find(~holesImage);
        if ~isempty(fidx)
            holesImage = ~imfill(holesImage,fidx(1),8);
            hR = regionprops(logical(holesImage),'PixelIdxList');
            holeStack = [];

            for h = 1:numel(hR)
                tmpH = zeros(size(holesImage));
                tmpH(hR(h).PixelIdxList) = 1;
                tmpH = imclose(tmpH,strel('disk',5,0));
                holeStack(:,:,h) = tmpH;
            end


        end
        if isempty(holeStack)
            holesImage = zeros(size(holesImage));
        else
            holesImage = any(holeStack,3);
        end
       




        %%%%%%%%%%%%%%%%%%%%%%%%%
        % make the color image
        %cR = regionprops(cLoop);
        cLoopLabel = bwlabel(cLoop);
        cLoopLabel = cLoopLabel + 1;
        resetIdx = find(useM == 0);
        cLoopLabel(resetIdx) = 0;
        rgbC = label2rgb(cLoopLabel);
        rgbC = double(rgbC)/255;
    catch ME
        getReport(ME);
    end
end