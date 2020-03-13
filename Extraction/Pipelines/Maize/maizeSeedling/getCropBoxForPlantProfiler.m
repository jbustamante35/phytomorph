function [qrCropBox,M] = getCropBoxForPlantProfiler(Lab,defaultAreaMinSize,internalBoxExpansionSize,reSZ,padValue,aThreshold,closeValue)
        % closeValue = 51;
        
        M = getMaskFromLab(Lab,aThreshold,closeValue);
        
        R = regionprops(M);
        box = R(1).BoundingBox;
        pad = round(padValue*reSZ);
        box(1:2) = box(1:2) - pad;
        box(3:4) = box(3:4) + 2*pad;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the boxes within the QR sheet
        fM = imfill(M,'holes');
        % remove small objects
        fM = bwareaopen(fM,defaultAreaMinSize*reSZ);
        % expand the internal boxes by 50 pixels
        fM = imdilate(fM,ones(round(internalBoxExpansionSize*reSZ)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the largest object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = regionprops(fM,'Area','BoundingBox');
        [J,midx] = max([R.Area]);
        if isempty(R(midx).BoundingBox)
            return
        end
        qrCropBox = R(midx).BoundingBox;
end
