function [result] = seedPlots(fileName,oPath,disp)
    try
    % inputs
        rimArea = [];
        result = 0;
    
        
        if nargin == 2;disp = false;end
        
        % set to plot edges and corner points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plotEdgeLines = false;
        plotCornerPoints = false;
        
        % default area filter size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        areaFilter = 30; % OLD before May 28, 2020
        areaFilter = 100; % NEW after May 28, 2020
        cardSquareSize = 1000;
        %cardSquareSize = 2244/2;
        
        % init table for phenotypes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        seedStatsTable = table();
        seedObjectTable = table();
        
        % config
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the file parts
        [~,nm,ext] = fileparts(fileName);
        
        fprintf(['start: read, color map and threshold data.\n']);tic
        % load data and init process
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image and convert color spaces
        I = imread(fileName);
        % convert RGB -> LAB
        LAB = rgb2lab(I);
        % get the b-channel
        blueYellow = LAB(:,:,3);
        % brightness channel
        brightness = bindVec(LAB(:,:,1));
        % brightness mask
        bMask = brightness < graythresh(brightness);
        % normalize the b-channel
        blueYellow = bindVec(blueYellow);
        M = blueYellow < graythresh(blueYellow);
        % get the seed mask from normalized channel
        seedMask = M;
        
        % fill holes on mask copy
        %blueMASK = imfill(M.*bMask,'holes');
        blueMASK = imfill(M,'holes');
        % THIS MIGHT NOT WORK FOR LAB PHOTOS
        blueMASK = imclearborder(blueMASK);

        blueMASK = bwlarge(blueMASK);
        fprintf(['end: read, color map and threshold data:' num2str(toc) ' \n']);


        fprintf(['start: find card edges.\n']);tic
        % find the card corners
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the edge
        E = edge(blueMASK);
        %E = edge(blueYellow);
        %E = bwopen(E,300);
        % perform hough transform
        %[H,T,R] = hough(E,'Theta',linspace(-90,89.99,1000));
        [H,T,R] = hough(E);
        % find peaks in hough
        %P  = houghpeaks(H,8,'threshold',ceil(0.3*max(H(:))),'NHoodSize',[201 51]);
        
        P  = houghpeaks(H,6,'threshold',ceil(0.1*max(H(:))));
        cnt = 1;
        for e = 1:size(P,1)
            ln = houghlines(E,T,R,P(e,:),'FillGap',1500,'MinLength',200);
            if ~isempty(ln)
                dis = [];
                for l = 1:numel(ln)
                    dis(l) = norm(ln(l).point1 - ln(l).point2);
                end
                [~,midx] = max(dis);
                lines(cnt) = ln(midx);
                lineMeasure(cnt) = H(P(e,1),P(e,2));
                cnt = cnt + 1;
            end
        end
                
        %{
        %P  = houghpeaks(H,8);
        GAP = 300:100:1800;
        LEN = 500:100:1800;
        for g = 1:numel(GAP)
            for l = 1:numel(LEN)
                % find lines
                lines = houghlines(E,T,R,P(1,:),'FillGap',GAP(g),'MinLength',LEN(l));
                
                numLines(g,l) = numel(lines);
                if numel(lines) == 4
                    mid = [];
                    for r = 1:numel(lines)
                        mid(r,:) = mean([lines(r).point1;lines(r).point2],1);
                        mid(r,:) = [lines(r).theta,lines(r).rho];
                    end
                    lineSpacing(g,l) = min(pdist(mid(:,2)));
                    
                    
                    % get intersections of lines
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    nidx = nchoosek(1:numel(lines),2);
                    COR = [];
                    for e = 1:size(nidx,1)
                        M1 = convertLineForm(lines(nidx(e,1)));
                        M2 = convertLineForm(lines(nidx(e,2)));
                        tmpP = intersect(M1,M2);
                        if ~isempty(tmpP);COR = [COR;tmpP];end
                    end
                    numCorners(g,l) = size(COR,1);
                    
                    if (numLines(g,l) == 4 && numCorners(g,l) == 4)
                        break;
                    end
                    
                    %if lineSpacing(g,l) > 40
                        close all
                        imshow(I,[]);hold on
                        range = [-4000 4000];
                        N = 10;
                        for l = 1:numel(lines)
                            M = convertLineForm(lines(l));
                            drawM(M,range,N);
                            if plotEdgeLines
                                xyl = [lines(l).point1;lines(l).point2];
                                plot(xyl(:,1),xyl(:,2),'m','LineWidth',2);
                            end
                            
                             for r = 1:numel(lines)
                                mid(r,:) = mean([lines(r).point1;lines(r).point2],1);
                             end
                             plot(mid(:,1),mid(:,2),'g*');
                        end
                        waitforbuttonpress
                    %end
                    
                else
                    lineSpacing(g,l) = 0;
                    numCorners(g,l) = 0;
                    
                end
            end
        end
        [gap,lines] = find(numLines == 4 & numCorners == 4);
        %[gap,lines] = find(numLines == 4 & lineSpacing > 400 & numCorners == 4);
        %[gap,lines] = find(numLines == 4);
        lines = houghlines(E,T,R,P,'FillGap',GAP(gap(1)),'MinLength',LEN(lines(1)));
        fprintf(['end: find card edges:' num2str(toc) ' \n']);
        %}
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % draw lines - potential
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp
            close all
            imshow(I,[]);hold on
            range = [-4000 4000];
            N = 10;
            for l = 1:numel(lines)
                M = convertLineForm(lines(l));
                drawM(M,range,N);
                if plotEdgeLines
                    xyl = [lines(l).point1;lines(l).point2];
                    plot(xyl(:,1),xyl(:,2),'m','LineWidth',2);
                end
            end
            drawnow
        end


        % generate a set of 4 lines to make a square
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        sidx = nchoosek(1:numel(lines),4);
        for s = 1:size(sidx,1)
            % select out the temp 4 lines
            tline = lines(sidx(s,:));
            % temp measure of each line
            tMeasure = lineMeasure(sidx(s,:));



            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % draw lines - potential
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if disp
                close all
                imshow(I,[]);hold on
                range = [-4000 4000];
                N = 10;
                for l = 1:numel(tline)
                    M = convertLineForm(tline(l));
                    drawM(M,range,N);
                    if plotEdgeLines
                        xyl = [tline(l).point1;tline(l).point2];
                        plot(xyl(:,1),xyl(:,2),'m','LineWidth',2);
                    end
                end
                drawnow
            end
            %}

            fprintf(['start: order points for ' num2str(s) ':' num2str(size(sidx,1)) ' .\n']);tic
            % draw points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            nidx = nchoosek(1:numel(tline),2);
            COR = [];
            ANG = [];
            for e = 1:size(nidx,1)
                M1 = convertLineForm(tline(nidx(e,1)));
                M2 = convertLineForm(tline(nidx(e,2)));
                [P,a] = intersect(M1,M2);
                if ~isempty(P);COR = [COR;P];ANG = [ANG;a];end
            end



            if size(COR,1) == 4
                % order points
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                dist = sum(COR.*COR,2);
                [~,midx] = min(dist);
                dist = abs(COR(midx,2) - COR(:,2));
                [~,tidx] = sort(dist);
                midx = [midx;tidx(2)];
                [~,tidx] = sort(COR(:,1));
                tidx = setdiff(tidx,midx,'stable');
                midx = [midx;tidx];
                fprintf(['end: order points:' num2str(toc) ' \n']);

                % measuredistance
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                orderedCOR = COR(midx,:);
                



                movingPoints = COR(midx,:);
                fixedPoints = [[-cardSquareSize -cardSquareSize];[cardSquareSize -cardSquareSize];[-cardSquareSize cardSquareSize];[cardSquareSize cardSquareSize]];

                % fit transformation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                tform = fitgeotrans(movingPoints,fixedPoints,'projective');
                itform = invert(tform);


                newA = 2000^2;
                A = itform.transformPointsForward(fixedPoints);
                dA = [[A(2,:) - A(1,:)];[A(3,:) - A(1,:)]];
                orgA = det(dA);
                ratioA(s) = abs(newA/orgA);


                di(1) = norm(orderedCOR(1,:) - orderedCOR(4,:));
                di(2) = norm(orderedCOR(1,:) - orderedCOR(4,:));
                ed(1) = norm(orderedCOR(1,:) - orderedCOR(2,:));
                ed(2) = norm(orderedCOR(2,:) - orderedCOR(4,:));
                ed(3) = norm(orderedCOR(4,:) - orderedCOR(3,:));
                ed(4) = norm(orderedCOR(1,:) - orderedCOR(3,:));
                ued = mean(ed);
                eed = (2*ued.^2).^.5;



                tMeasure = tMeasure(midx);
                tline = tline(midx);
                angleMeasure = abs(90*4 - sum(abs(ANG)));
                diagonalMeasure = eed - mean(di);

                if (ued > 100) && (ratioA(s) < 10) && (newA/orgA > 0)
                    squareMeasure(s) = angleMeasure + diagonalMeasure;
                else
                    squareMeasure(s) = inf;
                end

                
                %{
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % draw lines - potential
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disp
                    CL = {'r','b','g','k'};
                    close all
                    imshow(I,[]);hold on
                    range = [-4000 4000];
                    N = 10;
                    for l = 1:numel(tline)
                        M = convertLineForm(tline(l));
                        drawM(M,range,N);
                        if plotEdgeLines
                            xyl = [tline(l).point1;tline(l).point2];
                            plot(xyl(:,1),xyl(:,2),CL{l},'LineWidth',2);
                        end
                    end
                    drawnow
                end


                % plot points
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disp
                    for e = 1:numel(midx)
                        pt = COR(midx(e),:);
                        text(pt(1),pt(2),num2str(e))
                    end
                    drawnow
                    %waitforbuttonpress
                end
                %}
                


            else
                squareMeasure(s) = inf;
            end

        
        end


        [~,midx] = min(abs(squareMeasure));
        lines = lines(sidx(midx,:));
        selA = ratioA(midx);
        if  selA > 15
            throw(MException('h:l5','area error'));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % draw lines - potential
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp
            close all
            imshow(I,[]);hold on
            range = [-4000 4000];
            N = 10;
            for l = 1:numel(lines)
                M = convertLineForm(lines(l));
                drawM(M,range,N);
                if plotEdgeLines
                    xyl = [lines(l).point1;lines(l).point2];
                    plot(xyl(:,1),xyl(:,2),'m','LineWidth',2);
                end
            end
            drawnow
        end


        fprintf(['start: order points.\n']);tic
        % draw points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        nidx = nchoosek(1:numel(lines),2);
        COR = [];
        for e = 1:size(nidx,1)
            M1 = convertLineForm(lines(nidx(e,1)));
            M2 = convertLineForm(lines(nidx(e,2)));
            P = intersect(M1,M2);
            if ~isempty(P)
                COR = [COR;P];
                %{
                imshow(I,[]);
                xyl = [lines(nidx(e,1)).point1;lines(nidx(e,1)).point2];
                plot(xyl(:,1),xyl(:,2),'r','LineWidth',2);
                xyl = [lines(nidx(e,2)).point1;lines(nidx(e,2)).point2];
                plot(xyl(:,1),xyl(:,2),'k','LineWidth',2);
                plot(P(1),P(2),'g*')
                drawnow
                waitforbuttonpress
                %}
            end
        end

        % order points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        dist = sum(COR.*COR,2);
        [~,midx] = min(dist);
        dist = abs(COR(midx,2) - COR(:,2));
        [~,tidx] = sort(dist);
        midx = [midx;tidx(2)];
        [~,tidx] = sort(COR(:,1));
        tidx = setdiff(tidx,midx,'stable');
        midx = [midx;tidx];
        fprintf(['end: order points:' num2str(toc) ' \n']);

        % plot points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp
            if true
                for e = 1:numel(midx)
                    pt = COR(midx(e),:);
                    text(pt(1),pt(2),num2str(e))
                end
            end
            drawnow
        end

        % geo-warp images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pidx = midx;
        pidx(3:4) = flip(pidx(3:4),1);
        idealCardMask = poly2mask(COR(pidx,1),COR(pidx,2),size(I,1),size(I,2));
        
        % make square card
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        refR = imref2d(size(I));
        movingPoints = COR(midx,:);
        fixedPoints = [[-cardSquareSize -cardSquareSize];[cardSquareSize -cardSquareSize];[-cardSquareSize cardSquareSize];[cardSquareSize cardSquareSize]];

        % fit transformation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        tform = fitgeotrans(movingPoints,fixedPoints,'projective');
        itform = invert(tform);

        % map check
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        muX = mean(movingPoints);
        
        A = [tform.transformPointsForward(muX + [1 0]);
        tform.transformPointsForward(muX + [0 1])];
        B = [tform.transformPointsForward(muX + [0 0]);
        tform.transformPointsForward(muX + [0 0])];
        centerArea = det(A-B);


        %{
        if  centerArea < 0 || centerArea > 5
            
            centerArea
            rimArea
            selA
            throw(MException('h:l5','area error'));
        end
        %}

        % warp image card into blue square
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        [warpedI,tranR]  = imwarp(I,refR,tform);
        warpedCardMask = imwarp(idealCardMask,tform);
        warpedSeedMask = imwarp(seedMask,tform);
        warpedBrightness = imwarp(LAB(:,:,1),tform);
        warpedAchannel = imwarp(LAB(:,:,2),tform);
        warpedBchannel = imwarp(LAB(:,:,3),tform);
        toc
        
        % could use 1) pre warped or post warped threshold values
        % here - calc post warped seed mask on card only
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        fidx = find(warpedCardMask(:) > .5);
        
        subsetB = bindVec(warpedBchannel(fidx));
        maskV = subsetB > graythresh(subsetB);
        postWarpedSeedMask = zeros(size(warpedSeedMask));
        postWarpedSeedMask(fidx) = maskV;
        
        % normalize brightness
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        warpedBrightness = bindVec(warpedBrightness);
        
        % find bright areas
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        brightMask = warpedBrightness > graythresh(warpedBrightness);
        warpedSeedMask_mod = ~warpedSeedMask;
        warpedSeedMask_mod = imclearborder(warpedSeedMask_mod);
        
        % find bright areas that are seed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        warpedBrightSeed = (warpedSeedMask_mod & brightMask);
        
        % warp co-ordinates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [wy,wx] = ndgrid(1:size(warpedI,1),1:size(warpedI,2));
        uwx = imwarp(wx,tranR,itform,'OutputView',refR);
        uwy = imwarp(wy,tranR,itform,'OutputView',refR);

        %{
        % display seed overlay on warped image
        close all
        out = flattenMaskOverlay(warpedI,warpedBrightSeed>.5,.5,'r');
        imshow(out,[]);
        %}
        
        % image is ~seedMask
        warpedSeedMask = bwlarge(logical(warpedSeedMask));
        
        % use brightness to remove shawdows for outdoor images
        seedImage = warpedBrightSeed.*warpedCardMask;
        
        % IF NO EXTRA BRIGHTNESS FILTER FOR SHAWDOWS - INDOORS ONLY
        seedImage = ~warpedSeedMask;
        
        % USE POST WARPED - CARD ONLY -TEST INDOORS
        seedImage = postWarpedSeedMask;
        
        % area transform on objects
        seedImage = bwareaopen(seedImage,areaFilter);
        
        % clear border on objects
        seedImage = imclearborder(seedImage);
        
        % need to remove small "dust" objects
        %seedImage = bwareaopen(seedImage,100);

        % clear numbers - default on right
        %seedImage(:,2000:end) = 0;
        
        % get regionprops
        R = regionprops(seedImage,'Area','Perimeter','PixelIdxList','Centroid','Eccentricity');
        % count on area
        sidxArea = count([R.Area]);
        % count on perimeter
        sidxPerimeter = count([R.Perimeter]);
        % count on perimeter
        sidxEccentricity = count([R.Eccentricity]);


        kp = (sidxArea == sidxPerimeter) & (sidxArea == sidxEccentricity);
        kp = (sidxArea == sidxPerimeter);
        R = R(kp);
        sidx = count([R.Area]);




        %ri = round(.75*(mean([R(sidxArea==1).Area])/pi).^.5);

        % cluster colors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % unique number cluster counts
        UQ = unique(sidx);
        % get the jet colormap
        colorMap = jet;

        cx = linspace(0,1,size(colorMap,1));
        cxi = linspace(0,1,numel(UQ));
        y = interp1(cx,colorMap,cxi);


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        out = warpedI;
        CENL = [];
        unwarpedCenter = [];
        maskPT = zeros(size(warpedSeedMask));
        for u = 1:numel(UQ)
            % if the single cluster
            if UQ(u) == 1
                % 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %make a mask for the uth cluster
                newMask = zeros(size(warpedSeedMask));
                % find the uth clusters
                singleSeedIdx = find(sidx==UQ(u));
                % color them
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for e = 1:numel(singleSeedIdx)
                    newMask(R(singleSeedIdx(e)).PixelIdxList) = 1;
                end
                % overlay the uth cluster
                %out = flattenMaskOverlay(out,newMask>.9,.3,y(u,:));
                out = flattenMaskOverlay(out,newMask>.9,.5,'r');
                % sample the color space
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sidx = find(newMask>.9);
                seedColorData = [];
                for k = 1:3
                    tmp = warpedI(:,:,k);
                    seedColorData = [seedColorData tmp(sidx)];
                end
                [~,colorCom,colorU,colorE] = PCA_FIT_FULL(double(seedColorData)/255,1,false);

                % fill seed object table
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subR = R(singleSeedIdx);
                seedObjectTable = struct2table(subR);
                
                
                % fill the phenotypes table
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                seedStatsTable.averageArea = mean([R(singleSeedIdx(e)).Area]);
                seedStatsTable.averagePerimeter = mean([R(singleSeedIdx(e)).Perimeter]);
                seedStatsTable.averageRed = colorU(1);
                seedStatsTable.averageGreen = colorU(2);
                seedStatsTable.averageBlue = colorU(3);
                
                % isolate the centers
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % distance transform the uth cluster
                dT = bwdist(~newMask);
                % smooth the distance transform
                dT = imfilter(dT,fspecial('average',[3 3]),'replicate');
                CP = (imdilate(dT,strel('disk',11)) == dT).*newMask;
                maskPT = maskPT + CP;
                nR = regionprops(logical(CP),'Centroid');
                for c = 1:numel(nR)
                    CENL = [CENL;nR(c).Centroid];

                    tmpCenter = [];
                    ds = (uwx - nR(c).Centroid(1)).^2 +  (uwy - nR(c).Centroid(2)).^2;
                    [tmpCenter(:,1),tmpCenter(:,2)] = find(ds == min(ds(:)));
                    unwarpedCenter = [unwarpedCenter;tmpCenter(1,:)];
                end
            end

            % not sure what this is
            % i think it is trying to make a count
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %ALT(u) = numel(nR)/numel(singleSeedIdx);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        



        %{
        dT = bwdist(~seedImage);
        maskPTUSE = imdilate(maskPT,strel('disk',11,0));
        gidx = find(maskPTUSE);
        maskPTUSE(gidx) = dT(gidx) + 1;
        IM = imreconstruct(double(dT),maskPTUSE);
        skel = bwmorph(seedImage,'skeleton',inf);
        [s(:,1),s(:,2)] = find(skel);
        imshow(out,[]);hold on
        plot(s(:,2),s(:,1),'r.');
        rad = 100;
        [LEG,pos] = makeLegend(y,rad,[]);
        CENLD = tform.transformPointsForward(CENL);
        %}

        close all
        finalImage = imwarp(out,tranR,itform,'OutputView',refR);
        
        if disp
            imshow(finalImage,[]);
            hold on;
            plot(unwarpedCenter(:,2),unwarpedCenter(:,1),'k.');
        end
        

        if ~isempty(oPath)
            % write phenotype table
            writetable(seedStatsTable,[oPath nm '.csv']);
            % write color data
            csvwrite([oPath nm '_color.csv'],seedColorData);
            % write image overlay
            imwrite(finalImage,[oPath nm '.jpg']);
        end
        
        if disp
            drawnow
            pause(.1);
        end
        
        
        
        
        %phenotypeGraph = [];
        %[phenotypeGraph] = generatePhenotypeNode(phenotypeGraph,data,name,linker);
        
        %{
        hold on

        for l = 1:numel(lines)
            M = convertLineForm(lines(l));
            drawM(M,range,N);
            if plotEdgeLines
                xyl = [lines(l).point1;lines(l).point2];
                plot(xyl(:,1),xyl(:,2),'m','LineWidth',2);
            end
        end

        if plotCornerPoints
            for e = 1:numel(midx)
                pt = COR(midx(e),:);
                plot(pt(1),pt(2),'g*');
                text(pt(1),pt(2),num2str(e))
            end
        end


        if plotCornerPoints
            for e = 1:numel(midx)
                pt = COR(midx(e),:);
                text(pt(1),pt(2),num2str(e))
            end
        end
        %}

        %{
        imshow(LEG);
        for e = 1:size(pos,1)
            %cntLAB = [num2str(UQ(e)) '-' num2str(ALT(e))];
            cntLAB = [num2str(ALT(e))];
            text(2*rad,pos(e),cntLAB,'BackgroundColor','white');
        end
        %plot(CENLD(:,1),CENLD(:,2),'k.');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        if ~isempty(oPath)
            writeStatus(oPath,nm,1);
        end
    catch ME
        if ~isempty(oPath)
            writeStatus(oPath,nm,0);
        end
        getReport(ME)
        result = ME;
        result = rimArea;
    end

end

function [M] = convertLineForm(line)
    v = line.point2 - line.point1;
    v = v / norm(v);
    p = line.point1;
    M = [diag(v) p'];[0 0 1];
end

function [] = drawM(M,range,N,CL)
    if nargin == 3;CL = 'm';end
    L = linspace(range(1),range(2),N);
    pt = (M*[L;L;ones(size(L))])';
    plot(pt(:,1),pt(:,2),CL,'LineWidth',2);
end

function [a] = angle(M1,M2)
    d1 =diag(M1);
    d2 =diag(M2);
    a = acos(d1'*d2)*180/pi;
end

function [P,a] = intersect(M1,M2)
    a = angle(M1,M2);
    P = [];
    if a > 50
        options = optimset('MaxFunEvals',10000);
        f1 = @(a)(M1*[a*ones(2,1);1]);
        f2 = @(b)(M2*[b*ones(2,1);1]);
        func = @(x)sum(abs(f1(x(1)) - f2(x(2))));
        %x = fminunc(func,[1;1]);
        initValues = [[1;1],[-1;1],[1;-1],[-1;-1]];
        for e = 1:size(initValues,2)
            [x(:,e),fval(e)] = fminsearch(func,initValues(:,e),options);
        end
        [~,midx] = min(fval);
        x = x(:,midx);
        p1 = f1(x(1));
        p2 = f2(x(2));
        P = [p1';p2'];
        P = mean(P,1);
        %plot(p1(1),p1(2),'g*');
        %plot(p2(1),p2(2),'g*');
    end
end

function [LEG,pos] = makeLegend(cV,rad,pad)
    LEG = [];
    pos = rad;
    for e = 1:size(cV,1)
        tmpM = zeros(2*rad+1);
        tmpM((rad+1),(rad+1)) = 1;
        d = bwdist(tmpM);
        tmpM = (d < rad);
        kidx = find(tmpM==1);
        tmpL = ones(size(tmpM));
        disk = zeros([size(tmpL) 3]);
        for k = 1:3
            K = tmpL;
            K(kidx) = cV(e,k);
            disk(:,:,k) = K;
        end
        LEG = [LEG;disk];
        if e > 1;pos = [pos;pos(end)+2*rad];end
    end
    %imshow(LEG,[]);
end

%{
    FilePath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/';
    FileList = {};
    FileExt = {'JPG'};
    FileList = fdig(FilePath,FileList,FileExt,1);
    oPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/output2/';
    mmkdir(oPath);
    for e = 1:numel(FileList)
        seedPlots(FileList{e},oPath);
        %drawnow
        %waitforbuttonpress
    end
%}