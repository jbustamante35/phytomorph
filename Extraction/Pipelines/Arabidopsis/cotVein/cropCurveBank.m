function [dataPoint] = cropCurveBank(fileName,curveBank,oTemplate,disp)


    % get the affine transformation
    [d] = generateAffineImageDomain(oTemplate);
    % dilate the model mask to give some room
    modMask = imdilate(oTemplate > .8,strel('disk',5,0));


    [pth,nm,ext] = fileparts(fileName);
    I = double(imread(fileName))/255;
    I = rgb2gray(I);

    % for each curve in the bank
    for c = 1:size(curveBank,1)
        % init the cot model
        dataPoint(c) = cotModel(fileName,curveBank(c,:));

        % get the reference point
        trans0 = curveBank(c,1:5);
        trans = curveBank(c,6:10);

        % sample the image in the model frame
        fixedImageInMovingFrame = transformImage(oTemplate,I,d,trans,trans0);

        % threshold the image 
        th = graythresh(fixedImageInMovingFrame);
        fixedMaskInMovingFrame = fixedImageInMovingFrame > th;


        fixedMaskInMovingFrame = modMask .* fixedMaskInMovingFrame;
        fixedMaskInMovingFrame = logical(fixedMaskInMovingFrame > .8);
        %fixedMaskInMovingFrame = bwlarge(fixedMaskInMovingFrame);

        dataPoint(c).alignedImage = fixedImageInMovingFrame;
        dataPoint(c).alignedMask = fixedMaskInMovingFrame;



        if disp
            close all
            out = flattenMaskOverlay(bindVec(fixedImageInMovingFrame),fixedMaskInMovingFrame);
            %out = flattenMaskOverlay(out,oTemplate > .8,.3,'b');
            %out = flattenMaskOverlay(bindVec(fixedImageInMovingFrame),modMask>.8);
            imshow(out,[]);
            drawnow
            %pause(.4);
        end
    end
end