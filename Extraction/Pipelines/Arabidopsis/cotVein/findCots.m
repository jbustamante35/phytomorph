function [curBank] = findCots(cotSeachRegions,numberGroups,maxCots,templateData,fixedData,dispImage,optiPara)
    
    % set the template data
    moving = templateData.image;
    movingMask = templateData.mask;
    % set the fixed data
    fixed = fixedData.image;
    fixedMask = fixedData.mask;

    % set the opti parameters
    swarmSize = optiPara.swarmSize;
    bounds = optiPara.bounds;

    % get the cetner pointof the fixed image
    cp = hsize(fixedData.image);
    % get the region props for the search regions
    Rx = regionprops(logical(cotSeachRegions),'PixelIdxList');

    % init the data bank
    curBank = [];

    parfor loc = 1:min(numel(Rx),maxCots)
        tm = clock;
        % make mask for each search regions
        zidx = [];
        z = zeros(size(cotSeachRegions));
        z(Rx(loc).PixelIdxList) = 1;
        [zidx(:,1),zidx(:,2)] = find(z==1);
        [~,cidx] = kmeans(zidx,numberGroups);

        % display flags
        dispi = false;
        dispBest = false;
        dispf = false;
        

        xBest = [];
        fBest = [];
        x0Bank = [];

        for c = 1:size(cidx,1)
            x0Bank(c,:) = [0 0 0 cidx(c,1)-cp(1) cidx(c,2)-cp(2)];
            [xBest(c,:),fBest(c)] = myImageAlign(x0Bank(c,:),moving,movingMask,fixed,fixedMask,swarmSize,bounds,dispf,dispi,dispBest,dispImage);
            fprintf(['Done searching for cot:' num2str(loc) ': region #:' num2str(c) '.\n']);
        end

        % find the best for the image
        [~,midx] = min(fBest);
        curBank(loc,:) = [x0Bank(midx,:) xBest(midx,:)];
        etm = etime(clock,tm);
        fprintf(['Done searching for cot:' num2str(loc) ':' num2str(etm) '.\n']);
    end
end