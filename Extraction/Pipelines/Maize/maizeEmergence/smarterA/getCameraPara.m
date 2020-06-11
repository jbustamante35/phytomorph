function [tform] = getCameraPara(I)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Running: Get camera parameter(s) 2.0 \n']);
    % conversion
    squareSizeInMM = 25.4;
    % if image I is char - then read I
    if ischar(I);I = imread(I);end
    % get the checkerboard mask
    [MASK,boundingBox] = getCheckerBoardMask_ver2(I);
    % crop out the checkerboard from the image
    subI = imcrop(I,boundingBox);
    % built-in checker board point finding
    [imagePoints,boardSize] = detectCheckerboardPoints(subI);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract 5 points at a time from the point list
    % make 5 groups of 5 points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kidx = zeros(25,1);
    tmpImagePoints = imagePoints;
    for iter = 1:5
        q1 = min(tmpImagePoints(:,2));
        [JUNK,sidx] = sort(abs(tmpImagePoints(:,2) - q1),'ascend');
        kidx(sidx(1:5)) = iter;
        tmpImagePoints(sidx(1:5),2) = inf;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort up to down
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L = [];
    for g = 1:5
        L = [L ,kidx==g];
        v(g) = mean(imagePoints(kidx==g,2));
    end
    [~,sidx] = sort(v);
    L = L(:,sidx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort left to right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iS = [];
    for g = 1:5
        fidx = L(:,g)==1;
        d = imagePoints(fidx,:);
        [~,sidx] = sort(d(:,1));
        d = d(sidx,:);
        iS = [iS;d];
    end
    imagePoints = iS;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the center of the checker board points
    imageCenterPoint = [mean([boundingBox(:,1) boundingBox(:,1)+boundingBox(:,3)]) , mean([boundingBox(:,2) boundingBox(:,2)+boundingBox(:,4)])];
    % center the points on the image center
    imagePoints = bsxfun(@minus,imagePoints,[size(subI,2)/2 size(subI,1)/2]);
    % center the points on checker board "center"
    imagePoints = bsxfun(@plus,imagePoints,imageCenterPoint);
    % generate points
    worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort the generated points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort up to down
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kidx = kmeans(worldPoints(:,2),5);
    L = [];
    for g = 1:5
        L = [L ,kidx==g];
        v(g) = mean(worldPoints(kidx==g,2));
    end
    [~,sidx] = sort(v);
    L = L(:,sidx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort left to right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iS = [];
    for g = 1:5
        fidx = L(:,g)==1;
        d = worldPoints(fidx,:);
        [~,sidx] = sort(d(:,1));
        d = d(sidx,:);
        iS = [iS;d];
    end
    worldPoints = iS;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the work
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tform = fitgeotrans(worldPoints,imagePoints,'projective');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end