function [] = processFieldEarImage(fileName,oPath,areaThresh)
    % get file parts
    [pth,nm,ext] = fileparts(fileName);
    % make output image names
    border_oFile = [oPath nm '_border.csv'];
    central_oFile = [oPath nm '_central.csv'];
    image_oFile = [oPath nm '_return.jpg'];
    % read image
    I = imread(fileName);
    
    %{
    if rand(1) > .5
        I = imrotate(I,180);
    end
    %}
    % convert to LAB color space
    LAB = rgb2lab(I);
    % get the brightness signal
    sig0 = LAB(:,:,1);
    % get the yello/blue ratio
    sig = LAB(:,:,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the ears
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bind Vec
    sig = bindVec(sig);
    % threshold the signal
    th = graythresh(sig);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the white board
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig0 = bindVec(sig0);
    th0 = graythresh(sig0);
    whiteBoard = sig0 > th0;
    whiteBoard = imclose(whiteBoard,strel('disk',21,0));
    whiteBoard = imfill(whiteBoard,'holes');
    whiteBoard = bwlarge(whiteBoard);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % area filter the corn
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    MASK = sig > th & whiteBoard;
    MASK = bwareaopen(MASK,areaThresh);
    
    R = regionprops(MASK,'BoundingBox','MajorAxis','MinorAxis','Area','Centroid');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % area filter the corn
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find cos/sin wave
    gridSites = 2;
    fracDpi = .05;
    filterSize = 0;
    windowSize = round(1200*fracDpi):round(25*fracDpi):round(1600*fracDpi);
    CHUNK = 10;
    fftFilerValue = 0;
    for i = 1:numel(R)
        subMask = imcrop(MASK,R(i).BoundingBox);
        subMask = imfill(subMask,'holes');
        subI = imcrop(LAB(:,:,3),R(i).BoundingBox);
        subI = imrotate(subI,-90);
        subMask = imrotate(subMask,-90);
        T(i,:) = measureKernelOnEar(subI,subMask,filterSize,windowSize,gridSites,CHUNK,fftFilerValue);
    end
    T = nanmean(T,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:numel(T)
        box = R(i).BoundingBox;
        X = linspace(0,box(3),1000);
        %pos = [.5*(2*box(1)+box(3)) .5*(2*box(2)+box(4))];
        pos = [.5*(2*box(1)) .5*(2*box(2)+box(4))];
        Y = 10*cos(X*2*pi/T(i));
        plotSig{i} = [X'+pos(1) Y'+pos(2)];
    end
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter the centroid by the position
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    topV = [];
    for i = 1:numel(R)
        topV(i) = R(i).Centroid(2) > size(I,1)/2;
    end
    
    deltaTOP = [];
    deltaBOT = [];
    for i = 1:numel(R)
        if topV(i)
            deltaTOP = [deltaTOP R(i).Centroid(1)];
        else
            deltaBOT = [deltaBOT R(i).Centroid(1)];
        end
    end
    deltaTOP = sort(deltaTOP);
    deltaBOT = sort(deltaBOT);
    uT = mean(diff(deltaTOP));
    uB = mean(diff(deltaBOT));
    
    
    
    if uT > uB
        topBorder = true;
    else
        topBorder = false;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display mode for images - from version 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    out = flattenMaskOverlay(I,MASK);
    imshow(out,[]);
    hold on
    
    borderName = {};
    centralName = {};
    
    borderArea = [];
    centralArea = [];
    
    borderLength = [];
    centralLength = [];
    
    borderWidth = [];
    centralWidth = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over each
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numel(R)
        
        rectangle('Position',R(i).BoundingBox,'EdgeColor','g');
        L1(1,:) = R(i).BoundingBox(1:2);
        L1(2,:) = R(i).BoundingBox(1:2) + [R(i).BoundingBox(3) 0];
        L1(:,2) = L1(:,2) - 40;
        
        L2(1,:) = R(i).BoundingBox(1:2);
        L2(2,:) = R(i).BoundingBox(1:2) + [0 R(i).BoundingBox(4)];
        L2(:,1) = L2(:,1) - 40;
        
        plot(L1(:,1),L1(:,2),'b') 
        plot(L2(:,1),L2(:,2),'b') 
        
        T1 = L1(1,:) - 50;
        text(T1(1),T1(2),num2str(round(R(i).MinorAxisLength)));
        
        T2 = L2(1,:) - 20;
        h = text(T2(1),T2(2),num2str(round(R(i).MajorAxisLength)));
        set(h,'Rotation',90+180);
        
        
        plot(plotSig{i}(:,1),plotSig{i}(:,2),'b')
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display code for first generation method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if topBorder     
            if topV(i)
                nCL = 'm';
                borderArea = [borderArea R(i).Area];
                borderLength = [borderLength R(i).MajorAxisLength];
                borderWidth = [borderWidth R(i).MinorAxisLength];
                borderName{end+1} = nm;
            else
                nCL = 'c';
                centralArea = [centralArea R(i).Area];
                centralLength = [centralLength R(i).MajorAxisLength];
                centralWidth = [centralWidth R(i).MinorAxisLength];
                centralName{end+1} = nm;
            end
        else
            if topV(i)
                nCL = 'c';
                centralArea = [centralArea R(i).Area];
                centralLength = [centralLength R(i).MajorAxisLength];
                centralWidth = [centralWidth R(i).MinorAxisLength];
                centralName{end+1} = nm;
            else
                nCL = 'm'; 
                borderArea = [borderArea R(i).Area];
                borderLength = [borderLength R(i).MajorAxisLength];
                borderWidth = [borderWidth R(i).MinorAxisLength];
                borderName{end+1} = nm;
            end
        end
        dX = [100 -100];
        text(T1(1)+dX(1),T1(2)+dX(2),num2str(i),'BackgroundColor',nCL);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
    end
    hold off
    saveas(gca,image_oFile);
     
     
    toSpool = false;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if toSpool
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        borderTable = table(borderName',borderArea',borderLength',borderWidth','VariableNames',{'Name','Area','Length','Width'});
        centralTable = table(centralName',centralArea',centralLength',centralWidth','VariableNames',{'Name','Area','Length','Width'});
        writetable(borderTable,border_oFile);
        writetable(centralTable,central_oFile);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    drawnow
end