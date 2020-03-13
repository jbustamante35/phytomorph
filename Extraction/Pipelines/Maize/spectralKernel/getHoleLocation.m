function [holeLocations,R_seeds] = getHoleLocation(spectralImage,disp)

    uI = mean(spectralImage,3);
    V = bindVec(mean(mean(spectralImage,3),2));
    th = graythresh(V);
    th = .8;
    W = V > th;
    W = repmat(W,[1 size(spectralImage,2)]);
    
    W = imdilate(W,strel('square',41));
    W2 = imdilate(W,strel('square',41));
    R_seeds = regionprops(~W2,'BoundingBox');
    %%
    close all
    R_div = regionprops(W,'boundingBox');
    R_corn = regionprops(W,'boundingBox');
    holeLocations = [];
    for e = 1:numel(R_div)
        tmp = imcrop(uI,R_div(e).BoundingBox);
        tmp = bindVec(tmp);
        holes = tmp < graythresh(tmp);
        holes = imclearborder(holes);
        holes = bwlarge(holes,2);
        holesR = regionprops(holes,'Centroid');
        for h = 1:numel(holesR)
            holeLocations = [holeLocations;[holesR(h).Centroid + R_div(e).BoundingBox(1:2)]]
        end
        out = flattenMaskOverlay(double(tmp),holes);
        if disp
            imshow(out,[]);
            drawnow
        end
    end

end