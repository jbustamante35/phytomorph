function [tasselM] = thresholdTasselImage(img, reSZ)
    %
    %
    %orig = img;
    n = 5;
    for i=1:n

        Lab = rgb2lab(img);
        isColor = mean(abs(Lab(:,:,2:3)),3);
        raw = imcomplement(isColor);

        % Normalize tassel image and make pvc mask
        pvcMask = makePVCMask(img);

        % Create tassel mask
        tasselM = makeTasselMaskBlur(raw, reSZ, 500);

        %tasselM = makeTasselMaskKmeans(img);

        % Fill background of img at pvcMask with the mean of ~tasselMask
        img = fillBackground(tasselM, pvcMask, img);

        fprintf(['done with loop:' num2str(i) ':' num2str(n) '\n']);
    end
    
    tasselM = makeTasselMaskKmeans(img);

end
   