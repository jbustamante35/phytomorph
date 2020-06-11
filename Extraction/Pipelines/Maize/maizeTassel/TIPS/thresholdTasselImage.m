function [tasselM] = thresholdTasselImage(img, reSZ)
    %
    %
    %orig = img;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = 5;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop for background removal
    for i=1:n
        
        % convert RGB->lab
        Lab = rgb2lab(img);
        
        % get the mean a,b channels 
        isColor = mean(abs(Lab(:,:,2:3)),3);
        
        % if mean is high
        raw = imcomplement(isColor);

        % Normalize tassel image and make pvc mask
        pvcMask = makePVCMask(double(img));

        % Create tassel mask
        tasselM = makeTasselMaskBlur(raw, reSZ, 500);

        %tasselM = makeTasselMaskKmeans(img);

        % Fill background of img at pvcMask with the mean of ~tasselMask
        img = fillBackground(tasselM, pvcMask, img);

        fprintf(['done with loop:' num2str(i) ':' num2str(n) '\n']);
    end
    
    tasselM = makeTasselMaskKmeans(img);

end
   