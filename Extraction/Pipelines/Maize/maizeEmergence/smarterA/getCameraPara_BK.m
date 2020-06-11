function [tform] = getCameraPara(I)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract 5 points at a time from the point list
    % make 5 groups of 5 points
    kidx = zeros(25,1);
    tmpImagePoints = imagePoints;
    for iter = 1:5
        q1 = min(tmpImagePoints(:,2));
        [JUNK,sidx] = sort(abs(tmpImagePoints(:,2) - q1),'ascend');
        kidx(sidx(1:5)) = iter;
        tmpImagePoints(sidx(1:5),2) = inf;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make a matrix of [25 by 5]
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L = [];
    for g = 1:5
        L = [L ,kidx==g];
        v(g) = mean(imagePoints(kidx==g,2));
    end
    [~,sidx] = sort(v);
    L = L(:,sidx);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort left to right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iS = [];
    for g = 1:5
        fidx = L(:,g)==1;
        d = imagePoints(fidx,:);
        [~,sidx] = sort(d(:,1));
        d = d(sidx,:);
        iS = [iS;d];
    end
    imagePoints = iS;
    %imagePoints = flipdim(sortrows(flipdim(imagePoints,2)),2);
    %imagePoints = sortrows(imagePoints);
    %imagePoints = sortrows(flipdim(imagePoints,2));
    
    imageCenterPoint = [mean([boundingBox(:,1) boundingBox(:,1)+boundingBox(:,3)]) , mean([boundingBox(:,2) boundingBox(:,2)+boundingBox(:,4)])];
    %imagePoints = flipdim(imagePoints,2);
    %{
    imagePoints = flipdim(imagePoints,2);
    nP = size(imagePoints,1).^.5;
    [n1 n2] = ndgrid(linspace(0,(nP-1)*squareSizeInMM,nP),linspace(0,(nP-1)*squareSizeInMM,nP));
    n1 = flipdim(n1,1);
    worldPoints = [n1(:) n2(:)];
    worldPoints = flipdim(worldPoints,2);
    %}
    imagePoints = bsxfun(@minus,imagePoints,[size(subI,2)/2 size(subI,1)/2]);
    imagePoints = bsxfun(@plus,imagePoints,imageCenterPoint);
    worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
    %worldPoints = flipdim(worldPoints,2);
     
    kidx = kmeans(worldPoints(:,2),5);
    L = [];
    for g = 1:5
        L = [L ,kidx==g];
        v(g) = mean(worldPoints(kidx==g,2));
    end
    [~,sidx] = sort(v);
    L = L(:,sidx);
    
    % my sort
    iS = [];
    for g = 1:5
        fidx = L(:,g)==1;
        d = worldPoints(fidx,:);
        [~,sidx] = sort(d(:,1));
        d = d(sidx,:);
        iS = [iS;d];
    end
    worldPoints = iS;
    
    %worldPoints = flipdim(sortrows(flipdim(worldPoints,2)),2);
    %worldPoints = sortrows(worldPoints);
    %worldPoints = sortrows(flipdim(worldPoints,2));
    
    %worldPoints = bsxfun(@minus,worldPoints,mean(worldPoints,1));
    %worldPoints = bsxfun(@plus,worldPoints,imageCenterPoint);
    %worldPoints = flipdim(worldPoints,2);
    %worldCenterPoint = mean(worldPoints,1);
    %worldPoints = bsxfun(@minus,worldPoints,worldCenterPoint);
    %cameraParameters = estimateCameraParameters(cat(3,imagePoints,imagePoints),worldPoints);%,'ImageSize',(SZ(1:2)));
    tform = fitgeotrans(worldPoints,imagePoints,'projective');
    worldPoints
    imagePoints
end