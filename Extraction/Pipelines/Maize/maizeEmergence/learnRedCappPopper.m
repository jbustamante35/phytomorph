function [] = learnRedCappPopper(masterPath)

    try
        
        
        blackList = {'analysis-12','analysis-34'};
        
        
        csvKeyWord_Hand = 'fileName';
        csvKeyWord_Hand = 'Handscore';
        csvKeyWord_Green = 'Green';

        % display options for the mask
        viewMASK = false;

        % display options for the PCA
        viewAstack = false;
        viewMov = false;
        viewPLT = true;
        NUMtoView = 10;

        totalToSample = 5000;
        lowerDIMS = 10;
        imageScale = .5;
        imageScale2 = .75;
        diskDilateSize = 15;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scan for images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Scan for images.\n']);
        FilePath = masterPath;
        FileList = {};
        FileExt = {'tif'};
        FileList = fdig(FilePath,FileList,FileExt,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scan for images - ordered by folder
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Scan for images, order by folder.\n']);
        setFileList = {};
        [setFileList] = orderFrom_gdig(FileList,setFileList);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scan for csv
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Scan for csv.\n']);
        FilePath = masterPath;
        whole_csvFileList = {};
        FileExt = {'csv'};
        whole_csvFileList = fdig(FilePath,whole_csvFileList,FileExt,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % hand scores
        % scan for only csv files that have the key word in the name
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Scan for hand scores.\n']);
        kp = [];
        for e = 1:numel(whole_csvFileList)
            [p,n] = fileparts(whole_csvFileList{e});
            kp(e) = contains(n,csvKeyWord_Hand) & ~contains(n,'BK') & ~contains(n,'BAD');
        end
        hand_csvFileList = whole_csvFileList(find(kp));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % leaf scores
        % scan for only csv files that have the key word in the name
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Scan for leaf scores.\n']);
        kp = [];
        for e = 1:numel(whole_csvFileList)
            [p,n] = fileparts(whole_csvFileList{e});
            kp(e) = contains(n,csvKeyWord_Green) & ~contains(n,'BK') & ~contains(n,'BAD');
        end
        leaf_csvFileList = whole_csvFileList(find(kp));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a raw list only for set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for s = 1:numel(setFileList)
            kp = [];
            for e = 1:numel(setFileList{s})
                kp(e) = contains(setFileList{s}{e},'raw');
            end
            setFileList{s} = setFileList{s}(find(kp));
            % order the setFileList
            N = [];
            for e = 1:numel(setFileList{s})
                [p,n] = fileparts(setFileList{s}{e});
                nidx = strfind(n,'_');
                n = n(1:(nidx(1)-1));
                N(e) = str2num(n);
            end
            [~,sidx ] = sort(N);
            setFileList{s} = setFileList{s}(sidx);
            fprintf(['Done processing set:' num2str(s) ':' num2str(numel(setFileList)) '\n']);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a raw list only
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:numel(FileList)
            kp(e) = contains(FileList{e},'raw');
        end
        rFileList = FileList(find(kp));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute and sample N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sFileList = rFileList(randperm(numel(rFileList)));
        sFileList = sFileList(1:totalToSample);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read images to take the mean
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Start:Loading images for creating the mask for the center cap.\n'])
        I = [];
        errorRM = [];
        modN = 10;
        for e = 1:numel(sFileList)
            try
                tmp = imread(sFileList{e});
                oSZ = size(tmp);
                I(:,:,:,e) = imresize(tmp,imageScale);
                if (mod(e,nodN) == 0)
                    fprintf(['Done reading:' num2str(e) ':' num2str(totalToSample) '\n']);
                end
            catch ME
                errorRM = [errorRM e];
            end
        end
        % remove the error images
        I(:,:,:,errorRM) = [];
        fprintf(['End:Loading images for creating the mask for the center cap.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % determine if the images are the right size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SZ = [];
        fprintf(['Start:Scaning image for proper size.\n']);
        for s = 1:numel(setFileList)
            % fileparts
            [pth,nm,ext] = fileparts(setFileList{s}{1});
            % find the slash operator
            fidx = strfind(pth,'/');
            % get the cell number
            cellNumber = pth((fidx(end)+1):end);
            % get the analysis number and the folder name 
            analysisNumber = pth((fidx(end-3)+1):(fidx(end-2)-1));
            %analysisNumber = pth((fidx(end-1)+1):(fidx(end)-1));
            % create hash for set
            hIDX = ['H_' DataHash([analysisNumber '_' cellNumber])];
            % read image
            tmp = imread(setFileList{s}{1});
            SZ(s,:) = size(tmp);
            % store hash
            toUse.(hIDX) = all(SZ(s,:) == [301 301 3]);
            fprintf(['Done getting size for set:' num2str(s) ':' num2str(numel(setFileList)) '\n']);
        end
        fprintf(['End:Scaning image for proper size.\n']);
        %toUse = all(SZ(:,1) == 301 & SZ(:,2) == 301 & SZ(:,3) == 3,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the mask for the center cap
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Start: Finding the cap mask. \n'])
        close all
        % get mean
        U = mean(I,4);
        % find the center mask
        Ui = U;
        HSV = rgb2hsv(Ui/255);
        centerMask = HSV(:,:,2) > graythresh(HSV(:,:,2));
        centerMask = bwlarge(centerMask);
        centerMask = imdilate(centerMask,strel('disk',diskDilateSize));
        if viewMASK
            imshow(centerMask,[]);
            waitforbuttonpress
        end
        centerMask = imresize(centerMask,oSZ(1:2));
        centerMask = imresize(centerMask,imageScale2);
        % make the rings masks
        dt = bwdist(~(centerMask > .3));
        RAD = linspace(1,max(dt(:)),10);
        maskIDX = {};
        for r = 1:(numel(RAD)-1)
            msk = dt >= RAD(r) & dt <= RAD(r+1);
            maskIDX{r} = find(msk);
        end
        % get sample index for the mask
        kidx = find(centerMask > .3);
        fprintf(['End: Finding the cap mask. \n'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        funcSample = @(X)sampleRings(X,maskIDX);
        funcSample = @(X)sampleMaskOnly(X,kidx);
        
        
        
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample stratagy for X based on Y
        fairSampleSetX1 = {};
        fairSampleSetX2 = {};
        fairSampleSetX3 = {};
        numPerCellPreGerm = 2;
        numPerCellPostGerm = 3;
        for e = 1:numel(hand_csvFileList)
            
            % get the file parts for each csv file
            [p,n] = fileparts(hand_csvFileList{e});
            % find the slash operator
            fidx = strfind(p,'/');
            % get the cell number
            %cellNumber = p((fidx(end)+1):end);
            % read the Y measurements
            pulseData = csvread(hand_csvFileList{e});
            % get the analysis number and the folder name 
            analysisNumber = p((fidx(end-1)+1):(fidx(end)-1));
            % get the analysis number and the folder name 
            tmpPath = [masterPath analysisNumber filesep];
            
            
           
            Z = [];
            for cell = 1:numel(pulseData)
                try
                    fprintf(['starting cell:' num2str(cell) ':' num2str(numel(pulseData)) '\n'])
                    % create hash for set
                    tPath = [tmpPath 'output/imageStacks' filesep num2str(cell) filesep];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % get the list of images for this cell over time
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    FilePath = tPath;
                    tmpList = {};
                    FileExt = {'tif'};
                    tmpList = fdig(FilePath,tmpList,FileExt,1);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % get the raw images only
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    kp = [];
                    for i = 1:numel(tmpList)
                        kp(i) = contains(tmpList{i},'raw');
                    end
                    tmpList = tmpList(find(kp));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % order the tmpList
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    N = [];
                    for i = 1:numel(tmpList)
                        [p,n] = fileparts(tmpList{i});
                        nidx = strfind(n,'_');
                        n = n(1:(nidx(1)-1));
                        N(i) = str2num(n);
                    end
                    [~,sidx ] = sort(N);
                    tmpList = tmpList(sidx);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % make the popped vector
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    isPopped = pulseData(cell) ~= 0;
                    z = zeros(1,numel(tmpList));
                    if isPopped
                        z(pulseData(cell)) = 1;
                        z = cumsum(z);
                    end
                    Z(cell,:) = z;


                    r1 = [];r2 = [];
                    r1 = randi(sum(z==0),numPerCellPreGerm,1);
                    if isPopped
                        %r2 = randi(sum(z==1),numPerCellPostGerm,1);
                        r2 = find(z==1);
                        r2 = r2(1:min(numPerCellPostGerm,numel(r2)));
                        zidx = find(z==1);
                        r2 = zidx(1) + r2;
                    end

                    for r = 1:numel(r1)
                        fairSampleSetX1{end+1} = tmpList{r1(r)};
                    end


                    for r = 1:numel(r2)
                        fairSampleSetX2{end+1} = tmpList{r2(r)};
                    end
                    
                    fairSampleSetX3{end+1} = {{tmpList{r1}},{tmpList{r2}}};
                  
                    
                    fprintf(['ending cell:' num2str(cell) ':' num2str(numel(pulseData)) '\n'])
                catch ME
                    getReport(ME);
                end
                
            end
            
            e
            
        end
        
        %}
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %{
        
        N = 3000;
        fairSampleSetX1 = fairSampleSetX1(randperm(numel(fairSampleSetX1)));
        fairSampleSetX2 = fairSampleSetX2(randperm(numel(fairSampleSetX2)));
        for b = 1%:numel(blackList)
            fairSampleSetX1(contains(fairSampleSetX1,blackList{b})) = [];
            fairSampleSetX2(contains(fairSampleSetX2,blackList{b})) = [];
        end
        N = min([N numel(fairSampleSetX1) numel(fairSampleSetX2)]);
        sub1 = fairSampleSetX1(1:N);
        sub2 = fairSampleSetX2(1:N);
        
        %{
        sFileList = {sub1{:} sub2{:}};
        [scores1,error1,errorL1,U1,E1,L1] = sampleAndDecompose(sub1,imageScale2,funcSample,lowerDIMS);
        [scores2,error2,errorL2,U2,E2,L2] = sampleAndDecompose(sub2,imageScale2,funcSample,lowerDIMS);
        %}
        [X1] = sampleImageList(sub1,imageScale2,funcSample);
        [X2] = sampleImageList(sub2,imageScale2,funcSample);
         
        U1 = mean(X1,2);
        
        
        % create color function to project the data into a colored image
        tmp = imread(sub1{1});
        oSZ = size(tmp);
        tmp = imresize(tmp,imageScale2);
        szI = size(tmp);
        colorImage = @(colorData)putDataIntoImageAtP(colorData,kidx,szI);
        
        
        % place the sampled data into the image with the 
        U1i = colorImage(U1);
        U1i = imresize(U1i,[301 301]);
        kUi = reshape(U1i,[prod(oSZ(1:2)) oSZ(3)]);
        kidx2 = kmeans(kUi,3);
        kidx2 = reshape(kidx2,[301 301]);
        
        ringMask = bwlarge(kidx2==3);
        capMask = bwlarge(kidx2==2);
        
        
        alignMask = imresize(centerMask,oSZ(1:2));
        alignMask = bwlarge(alignMask==1);
        aR = regionprops(alignMask);
        [a1,a2] = ndgrid(-aR.BoundingBox(3)/2:aR.BoundingBox(3)/2,-aR.BoundingBox(4)/2:aR.BoundingBox(4)/2);
        x.domain = [a1(:) a2(:) ones(size(a1(:)))];
        x.sz = size(a1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make model to align to
        fT = @(p,globData)[[eye(2) p'];[0 0 1]];
        cp = csize(U1i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fill in with nan
        for k = 1:size(U1i,3)
            tmp = U1i(:,:,k);
            tmp(tmp(:)==0) = NaN;
            U1i(:,:,k) = tmp;
        end
        modelP = funcQ(cp(1:2),fT,U1i);
        modelF = modelP.getF(x);
        
        
        maskP = funcQ(cp(1:2),fT,double(alignMask));
        maskF = maskP.getF(x);
        
        for tr = 1:size(X1,2)

            trialI = double(imread(sub1{tr}));
            
            %trialI = colorImage(X1(:,tr));
            %trialI = imresize(trialI,oSZ(1:2));
            trialP = funcQ(cp(1:2),fT,trialI);
            trialF = trialP.getF(x);

            res = modelP.morph(trialP,x,true);
            resF = res.getF(x);
            
            trialF_masked = bsxfun(@times,trialF,maskF);
            resF_masked = bsxfun(@times,resF,maskF);
            
            figure
            imshow([modelF trialF_masked resF_masked]/255,[]);
            drawnow
        end
        
        morph(obj,target,x)
        
        dX_range = [[50;0],[50;0]];
        dR_range = [0,0];
        dS_range = [[2;.5],[2;.5]];
        funcT = fwdT.makeTF(dX_range,dR_range,dS_range);
        
        
        
        
        
        tmp = imread(sub1{1});
        oSZ = size(tmp);
        tmp = imresize(tmp,imageScale2);
        szI = size(tmp);
        colorImage = @(colorData)putDataIntoImageAtP(colorData,kidx,szI);
        
        
        for e = 1:size(X1,2)
            tmp = colorImage(X1(:,e));
            imshow(tmp/255,[]);
            %title(num2str(sub2{e}))
            title(num2str(e))
            waitforbuttonpress
        end
        
        
        
        PCA_FIT_FULL_T2ws(X1,X2,[2 7 3],colorImage);
        %{
        scores = PCA_REPROJ_T(I,E,U);
        sims = PCA_BKPROJ_T(scores,E,U);
        %}
        
        
        
        
        %}
        



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read and sample the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Start:Loading images for creating the PCA.\n']);
        errorRM = [];
        % load first frame and pre-allocate data block
        tmp = imread(sFileList{1});
        oSZ = size(tmp);
        tmp = imresize(tmp,imageScale2);
        [tmpColorSample] = funcSample(tmp);
        I = zeros(numel(tmpColorSample),numel(sFileList));
        parfor e = 1:numel(sFileList)
            try
                tmp = imread(sFileList{e});
                oSZ = size(tmp);
                tmp = imresize(tmp,imageScale2);
                [tmpColorSample] = funcSample(tmp);
                I(:,e) = tmpColorSample;
                fprintf(['Done reading:' num2str(e) ':' num2str(totalToSample) '\n']);
                rm(e) = false;
            catch ME
                rm(e) = true;
                getReport(ME)
            end
        end
        %errorRM = find(all(I == 0),1);
        errorRM = find(rm);
        % remove the error images
        I(:,errorRM) = [];
        fprintf(['End:Loading images for creating the PCA.\n'])
        % decompose the samped image data
        [U,E,L] = PCA_FIT_FULL_Tws(I,lowerDIMS);
        cc = PCA_REPROJ_T(I,E,U);
        simc = PCA_BKPROJ_T(cc,E,U);
        errorL = mean(sum((I - simc).*(I - simc),1),2);
        errorS = std(sum((I - simc).*(I - simc),1),1,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %{
        [f,xi] = ksdensity(sum((I - simc).*(I - simc),1).^.5);
        setSelect = 1;
        imageScale = .75;
        tmpE = [];
        cnt = 1;
        skip = 5;
        close all
        for fr = 1:10:numel(setFileList{setSelect})
            fileName = setFileList{setSelect}{fr};
            [~,tmpE(cnt)] = gatherPCAdata(fileName,imageScale,funcSample,U,E);
            tmpE(cnt) = errorP(tmpE(cnt)^.5);
            cnt = cnt + 1;
            plot(tmpE)
            drawnow
            fprintf(['done with ' num2str(fr) ':' num2str(numel(setFileList{setSelect})) '\n']);
        end
        %}
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute and sample N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sFileList = rFileList(randperm(numel(rFileList)));
        sFileList = sFileList(1:totalToSample);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack viewing for display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if viewAstack
            % permute the sets
            viewIDX = randperm(numel(setFileList));
            viewIDX = viewIDX(1:NUMtoView);
            viewIDX = 13;
            scores = [];
            % for each index
            for s = 1:numel(viewIDX)
                try
                    % set the selection indx
                    setSelect = viewIDX(s);
                    scores = [];
                    % get the index information
                    [p,n] = fileparts(setFileList{setSelect}{1});
                    fidx = strfind(p,'/');
                    cellNumber = p((fidx(end)+1):end);
                    analysisNumber = p((fidx(end-3)+1):(fidx(end-2)-1));
                    sidx = [];
                    for e = 1:numel(csvFileList)
                        sidx(e) = contains(csvFileList{e},analysisNumber);
                    end
                    sidx = find(sidx);
                    pulseData = csvread(csvFileList{sidx});
                    pulseData = pulseData(str2num(cellNumber));
                    pulseData_T = zeros(1,numel(setFileList{setSelect}));
                    if pulseData ~= 0
                        pulseData_T(pulseData) = 1;
                    end

                    scores = [];
                    for e = 1:numel(setFileList{setSelect})
                        
                        tmp = imread(setFileList{setSelect}{e});
                        tmp = imresize(tmp,imageScale2);
                        [tmpColorSample] = funcSample(tmp);
                        
                        %{
                        scores(:,e) = PCA_REPROJ_T(tmpColorSample,E,U);
                        tmpErr = tmpColorSample - PCA_BKPROJ_T(scores(:,e),E,U);
                        error(e) = norm(tmpErr);
                        %}
                        
                        %scores(:,e) = scores(:,e).*L.^-.5;
                        if viewMov
                            imshow(tmp,[]);
                            title(num2str(e))
                            drawnow
                        end
                    end

                    if viewPLT
                        plot(scores')
                        hold all
                        plot(pulseData_T,'k')
                        drawnow
                        waitforbuttonpress
                        hold off
                    end
                catch
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack the data for training
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toLoadX = true;
        toLoadY = true;
        toLoadCellPaths = true;
        
        
        if toLoadX;X = {};X1={};X2={};end
        if toLoadY;Y1 = [];end
        if toLoadY;Y2 = {};end
        if toLoadCellPaths;CP = [];end
        
        cnt = 1;
        viewLoad = false;
        % loop over each csv file
        for l = 1:numel(hand_csvFileList)
            fprintf(['Starting csv list:' num2str(l) ':' num2str(numel(hand_csvFileList)) '\n'])
            
            % get the file parts for each csv file
            [p,n] = fileparts(hand_csvFileList{l});
            % find the slash operator
            fidx = strfind(p,'/');
            % get the cell number
            %cellNumber = p((fidx(end)+1):end);
            % read the Y measurements
            pulseData = csvread(hand_csvFileList{l});
            % get the analysis number and the folder name 
            analysisNumber = p((fidx(end-1)+1):(fidx(end)-1));
            % get the analysis number and the folder name 
            %analysisNumber = pth((fidx(end-4)+1):(fidx(end-2)-1));
            % create hash for set
            hIDX = ['H_' DataHash([analysisNumber '_' cellNumber])];
            
            % tmp path for the 
            tmpPath = [masterPath analysisNumber filesep];
            CP(l).highPath = tmpPath;
            
            % for each cell
            CELLpathList = {};
            tmpScore = {};
            tmpY2 = {};
            tmpY1 = [];
            tmpLP = {};
            for cell = 1:numel(pulseData)
                fprintf(['Starting cell list:' num2str(cell) ':' num2str(numel(pulseData)) '\n'])
                try
                    % get the hash data
                    hIDX = ['H_' DataHash([analysisNumber '_' num2str(cell)])];

                    if toUse.(hIDX)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % construct the path to the image stacks
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        tPath = [tmpPath 'output/imageStacks' filesep num2str(cell) filesep];
                        CELLpathList{cell} = tPath;
                        
                        if toLoadX | toLoadY
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % get the list of images for this cell over time
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            FilePath = tPath;
                            tmpList = {};
                            FileExt = {'tif'};
                            tmpList = fdig(FilePath,tmpList,FileExt,1);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % get the raw images only
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            kp = [];
                            for e = 1:numel(tmpList)
                                kp(e) = contains(tmpList{e},'raw');
                            end
                            tmpList = tmpList(find(kp));
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % order the tmpList
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            N = [];
                            for e = 1:numel(tmpList)
                                [p,n] = fileparts(tmpList{e});
                                nidx = strfind(n,'_');
                                n = n(1:(nidx(1)-1));
                                N(e) = str2num(n);
                            end
                            [~,sidx ] = sort(N);
                            tmpList = tmpList(sidx);



                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % make the popped vector
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            isPopped = pulseData(cell) ~= 0;
                            z = zeros(1,numel(tmpList));
                            if isPopped
                                z(pulseData(cell)) = 1;
                                z = cumsum(z);
                            end
                        end
                        
                        
                        
                        if toLoadX
                            scores = [];
                            scoresPOST = [];
                            scoresPRE = [];
                            scores2 = [];
                            error = [];
                            errorPRE = [];
                            errorPOST = [];
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % read the image stack data and project to the scores
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            parfor tm = 1:numel(tmpList)
                                %fprintf(['Done with image:' num2str(tm) ':' num2str(numel(tmpList)) '@cell:' num2str(cell) ':' num2str(numel(pulseData)) '\n']);
                                % read the image
                                tmp = imread(tmpList{tm});
                                % resize the image
                                tmp = imresize(tmp,imageScale2);
                                % convert to lab for second signal extraction
                                tmpLAB = rgb2lab(tmp);
                                % reshape the tmpLAB
                                tmpSZ = size(tmpLAB);
                                tmpLAB = reshape(tmpLAB,[prod(tmpSZ(1:2)) tmpSZ(3)]);
                                sig2 = mean(tmpLAB(centerMask==1,:),1);
                                if viewLoad
                                    imshow(double(tmp).*centerMask/255,[]);
                                    title(num2str(tm))
                                    drawnow
                                end
                                [tmpColorSample] = funcSample(tmp);
                                
                                % for single vector space
                                scores(:,tm) = PCA_REPROJ_T(tmpColorSample,E,U,false);
                                tmpErr = tmpColorSample - PCA_BKPROJ_T(scores(:,tm),E,U);
                                error(tm) = norm(tmpErr);
                                
                                %{
                                % for double vector space(s)
                                scoresPRE(:,tm) = PCA_REPROJ_T(tmpColorSample,E1,U1,false);
                                scoresPOST(:,tm) = PCA_REPROJ_T(tmpColorSample,E2,U2,false);
                                tmpErrPRE = tmpColorSample - PCA_BKPROJ_T(scoresPRE(:,tm),E1,U1);
                                tmpErrPOST = tmpColorSample - PCA_BKPROJ_T(scoresPOST(:,tm),E2,U2);
                               
                                 
                                errorPRE(tm) = norm(tmpErrPRE);
                                errorPOST(tm) = norm(tmpErrPOST);
                                %}
                                
                                scores2(:,tm) = sig2';
                            end
                            % store the scores
                            %X{cnt} = [scores;error];
                            % for single vectors
                            tmpScore{cell} = scores;
                            % for double vectors
                            %tmpScore{cell} = [scoresPRE;scoresPOST];
                            
                            tmpScore2{cell} = scores2;
                            
                            % for single error
                            tmpScore3{cell} = error;
                            
                            % for double error
                            %tmpScore3{cell} = [errorPRE;errorPOST];
                        end
                        
                        
                        if toLoadY
                            % store the isPopped vector
                            tmpY1(cell) = categorical(isPopped);
                            % store the time-popped vector
                            tmpY2{cell} = categorical(z);
                            % image stack information
                            tmpLP{cell} = tPath;
                        end
                        
                        
                    
                    end
                catch ME
                    ME
                end
            end
            
            
            for cell = 1:numel(tmpScore)
                X{cnt} = tmpScore{cell};
                X2{cnt} = tmpScore2{cell};
                X3{cnt} = tmpScore3{cell};
                Y2{cnt} = tmpY2{cell};
                Y1 = [Y1;tmpY1'];
                LP{cnt} = tmpLP{cell};
                cnt = cnt + 1;
            end
            
            
            CP(l).pathList = CELLpathList;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        blackList = 'analysis-12';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % process signals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sz = 11;
        cntT = 1;
        nX = {};
        nY = {};
        nLP = {};
        nY_cat = [];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create new X and new Y
        for e = 1:numel(Y2)
            if ~isempty(LP{e})
                % do not add the black listed data
                if ~contains(LP{e},blackList)

                    
                    osig = X{e};
                    esig = X3{e};
                    %osig = osig(sigsToUse,:);
                    
                    % error for double
                    %esig = esig(1:2,:);
                    
                    
                    
                    % cat error onto osig
                    %osig = [osig ; esig];
                    
                    %try zscore normalize
                    %sig = zscore(osig,1,2);
                    
                    % try filter
                    %osig = imfilter(osig,ones(1,sz)/sz,'replicate');
                    
                    % set the first frame to zero
                    sig = bsxfun(@minus,osig,osig(:,1));
                    
                    
                    %[L;errorL].^-.5
                    
                    %%%%%%%%%%% NOTE OUT
                    % eigenvalue normalize - with error
                    %sig = bsxfun(@times,sig,);
                    
                    % eignevalue normalize - without error
                    %sig = bsxfun(@times,sig,[L].^-.5);
                    
                    % normalize by std for double
                    %sig = bsxfun(@times,sig,[L1;L2;errorL1;errorL2].^-.5);
                    
                    % look at other work ups like rates - etc
                    %[J,sidx] = sort(sum(sig,2),'descend');
                    %sig = sig(sidx,:);
                    %sig = [sum(sig,1);std(sig,1,1)];
                    %nX{e} = [nX{e} ; gradient(sig)];

                    nX{cntT} = sig;
                    nY{cntT} = Y2{e};
                    nY_cat(cntT) = Y1(e);
                    nLP{cntT} = LP{e};
                    cntT = cntT + 1;
                    %nX{e} = [nX{e} ; X{e}];
                end
            end
        end
        
        %{
        for e = 1:numel(Y2)
            nX{e} = convertToStrain(X{e},11,10);
        end
        %}
        %{
        for e = 1:numel(Y2)
            nX{e} = bsxfun(@minus,nX{e},min(nX{e},[],2));
            nX{e} = bsxfun(@times,nX{e},max(nX{e},[],2).^-1);
        end
        %}
        
        
        % create new Y2
        for e = 1:numel(Y2)
            fidx = find(double(Y2{e}) == 2);
            z = zeros(size(Y2{e}));
            if ~isempty(fidx)
                fidx = fidx(1);
                z(fidx) = 1;
                z = imdilate(z,strel('disk',3,0));
                fidx = find(z==1);
                z(fidx(end)+1:end) = 2;
            end
            nY2{e} = categorical(z);
        end
        
        
        
        kp = [];
        for e = 1:numel(nX)
            if (size(nX{e},2) == size(nY{e},2))
                %if ~all(double(nY{e})==1)
                    kp = [kp e];
                %end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        loopScan = false;
        if loopScan
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % training scanner over pca, error and internal states for
            % lstm-network
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            W = 15;
            D = 10;

            metricStatsCourse = {};
            metricStatsFine = {};
            trainMag = 5;
            % for each error
            for useError = 0:1
                % for numberof PCs
                for numPC = 1:10
                    % for internal states
                    for states = 1:10
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % stack the data for this round of testing
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        useX = {};
                        for e = 1:numel(nX)
                            tmpX = nX{e}(1:numPC,:);
                            if useError;tmpX = [tmpX;nX{e}(end,:)];end
                            useX{e} = tmpX;
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % create cross validation partition
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        cv = cvpartition(numel(kp),'HoldOut',.2);
                        trainIDX = find(cv.training);
                        testIDX = find(cv.test);
                        % get the training set
                        trainX = useX(kp(trainIDX));
                        trainY = nY(kp(trainIDX));
                        % get the test set
                        testX = useX(kp(testIDX));
                        testY = nY(kp(testIDX));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make fine net data - training
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fineXTrain = {};
                        fineYTrain = {};
                        for tr = 1:numel(trainIDX)
                            fidx = find((double(trainY{tr})-1) == 1);
                            centerP = 0;
                            if ~isempty(fidx)
                                dD = randi(2*D+1,1)-(1+D);
                                centerP = fidx(1) + dD;
                            end
                            widx = (centerP-W):(centerP+W);
                            if (widx(1) > 1) && (widx(end) < numel(trainY{tr}))
                                cutX = trainX{tr}(:,widx);
                                cutY = trainY{tr}(widx);
                                fineXTrain{end+1} = cutX;
                                fineYTrain{end+1} = cutY;
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make fine net data - testing
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        fineXTest = {};
                        fineYTest = {};
                        for tr = 1:numel(testIDX)
                            fidx = find((double(testY{tr})-1) == 1);
                            centerP = 0;
                            if ~isempty(fidx)
                                dD = randi(2*D+1,1)-(1+D);
                                centerP = fidx(1) + dD;
                            end
                            widx = (centerP-W):(centerP+W);
                            if (widx(1) > 1) && (widx(end) < numel(testY{tr}))
                                cutX = testX{tr}(:,widx);
                                cutY = testY{tr}(widx);
                                fineXTest{end+1} = cutX;
                                fineYTest{end+1} = cutY;
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % training for fine net
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        options = trainingOptions('sgdm',...
                                'InitialLearnRate',.001,...
                                'MaxEpochs',10*trainMag,...
                                'MiniBatchSize',128,...
                                'Shuffle','every-epoch',...
                                'Verbose',true,...
                                'ExecutionEnvironment','cpu',...
                                'Plots','training-progress');
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make layers
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        whenPoppedNetlayers = [ ...
                            sequenceInputLayer(size(trainX{1},1))
                            lstmLayer(states)
                            fullyConnectedLayer(2)
                            softmaxLayer
                            classificationLayer];
                        [whenPoppedNET,ti] = trainNetwork(fineXTrain,fineYTrain,whenPoppedNetlayers,options);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % grade the fine
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        met = [];
                        preY = whenPoppedNET.predict(fineXTest);
                        for te = 1:numel(preY)
                            sig1 = (double(fineYTest{te}) - 1);
                            sig2 = double(preY{te}(2,:));
                            f1 = find(sig1 > .5);
                            f2 = find(sig2 > .5);
                            if isempty(f1);f1 = 0;end
                            if isempty(f2);f2 = 0;end
                            met(te,1) = norm(sig1 - sig2);
                            met(te,2) = abs(f1(1) - f2(1));
                        end
                        metricStatsFine{states,useError+1,numPC} = met;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % training for course net
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        options = trainingOptions('sgdm',...
                                'InitialLearnRate',.0001,...
                                'MaxEpochs',10*trainMag,...
                                'MiniBatchSize',128,...
                                'Shuffle','every-epoch',...
                                'Verbose',true,...
                                'ExecutionEnvironment','cpu',...
                                'Plots','training-progress');
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make layers
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        whenPoppedNetlayers = [ ...
                            sequenceInputLayer(size(trainX{1},1))
                            lstmLayer(states)
                            fullyConnectedLayer(2)
                            softmaxLayer
                            classificationLayer];
                        [whenPoppedNET,ti] = trainNetwork(trainX,trainY,whenPoppedNetlayers,options);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % grade the course
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        met = [];
                        preY = whenPoppedNET.predict(testX);
                        for te = 1:numel(preY)
                            sig1 = (double(testY{te}) - 1);
                            sig2 = double(preY{te}(2,:));
                            f1 = find(sig1 > .5);
                            f2 = find(sig2 > .5);
                            if isempty(f1);f1 = 0;end
                            if isempty(f2);f2 = 0;end
                            met(te,1) = norm(sig1 - sig2);
                            met(te,2) = abs(f1(1) - f2(1));
                        end

                        metricStatsCourse{states,useError+1,numPC} = met;

                        close all force



                        coursePlot = [];
                        for s = 1:size(metricStatsCourse,1)
                            for err = 1:size(metricStatsCourse,2)
                                for pc = 1:size(metricStatsCourse,3)
                                    tmp = mean(metricStatsCourse{s,err,pc},1);
                                    if isempty(tmp);tmp = [0 0];end
                                    coursePlot(s,err,pc,:) = tmp;
                                end
                            end
                        end

                        finePlot = [];
                        for s = 1:size(metricStatsFine,1)
                            for err = 1:size(metricStatsFine,2)
                                for pc = 1:size(metricStatsFine,3)
                                    tmp = mean(metricStatsFine{s,err,pc},1);
                                    if isempty(tmp);tmp = [0 0];end
                                    finePlot(s,err,pc,:) = tmp;
                                end
                            end
                        end

                        figure;
                        hold all
                        for pc = 1:size(finePlot,3)
                            plot(squeeze(finePlot(:,1,pc,2)),'DisplayName',['pc' num2str(pc) '-f']);

                            plot(squeeze(coursePlot(:,1,pc,2)),'--','DisplayName',['pc' num2str(pc) '-c']);
                         end
                        legend('show');
                        xlabel('internal states');
                        ylabel('delta Frames');
                        hold on

                    end



                end

            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fine signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W = 15;
        D = 10;
        %%%%%%%%%%%%%%%%%%%%%%%%
        % signals to use
        sigsToUse = [[1:7],[10:17]]';
        % signals to use for single not double
        sigsToUse = [1:7]';
        useError = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare X and Y for fine net
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack the data for this round of testing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        useX = {};
        numPC = 7;
        for e = 1:numel(nX)
            tmpX = nX{e}(sigsToUse,:);
            if useError;tmpX = [tmpX;nX{e}(end,:)];end
            useX{e} = tmpX;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        useY = nY;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare X and Y for fine net      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make fine net data - testing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fineXTrain = {};
        fineYTrain = {};
        for tr = 1:numel(useX)
            fidx = find((double(useY{tr})-1) == 1);
            centerP = 0;
            if ~isempty(fidx)
                dD = randi(2*D+1,1)-(1+D);
                centerP = fidx(1) + dD;
            end
            widx = (centerP-W):(centerP+W);
            if (widx(1) > 1) && (widx(end) < numel(useY{tr}))
                cutX = useX{tr}(:,widx);
                cutY = useY{tr}(widx);
                fineXTrain{end+1} = cutX;
                fineYTrain{end+1} = cutY;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare X and Y for fine net
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % training for fine net
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        options = trainingOptions('sgdm',...
                'InitialLearnRate',.001,...
                'MaxEpochs',100,...
                'MiniBatchSize',128,...
                'Shuffle','every-epoch',...
                'Verbose',true,...
                'ExecutionEnvironment','cpu',...
                'Plots','training-progress');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        whenPoppedNetlayers = [ ...
            sequenceInputLayer(size(fineXTrain{1},1))
            lstmLayer(7)
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];
        [whenPoppedNET_fine] = trainNetwork(fineXTrain,fineYTrain,whenPoppedNetlayers,options);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fine signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% course signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare X and Y for course net
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % stack the data for this round of testing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        useX = {};
        numPC = 7;
        sigsToUse = [[1:7],[10:17]]';
        for e = 1:numel(nX)
            tmpX = nX{e}(sigsToUse,:);
            if useError;tmpX = [tmpX;nX{e}(end,:)];end
            useX{e} = tmpX;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        useY = nY;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % training for course net
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        options = trainingOptions('sgdm',...
                'InitialLearnRate',.001,...
                'MaxEpochs',10,...
                'MiniBatchSize',128,...
                'Shuffle','every-epoch',...
                'Verbose',true,...
                'ExecutionEnvironment','cpu',...
                'Plots','training-progress');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        whenPoppedNetlayers = [ ...
            sequenceInputLayer(size(useX{1},1))
            lstmLayer(7)
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];
        [whenPoppedNET_course] = trainNetwork(useX,useY,whenPoppedNetlayers,options);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        
        
        
        % save algo
        sPath = '/mnt/snapper/nate/redCapAlgo/';
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        options = trainingOptions('sgdm',...
                'InitialLearnRate',1,...
                'MaxEpochs',100,...
                'MiniBatchSize',128,...
                'Shuffle','every-epoch',...
                'Verbose',true,...
                'ExecutionEnvironment','cpu',...
                'Plots','training-progress');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        isPoppedNetlayers = [ ...
            sequenceInputLayer(size(useX{1},1))
            lstmLayer(7,'OutputMode','last')
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];
        for e = 1:numel(nY)
            catY(e) = categorical(nY{e}(end));
        end
        isPoppedNET = trainNetwork(useX,catY',isPoppedNetlayers,options);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        for analysis = 1:numel(CP)
            tmpOpath = [CP(analysis).highPath 'output' filesep];
            [greenPoppedFrame] = automateGreenHandScore(CP(analysis).pathList,imageScale2,centerMask,tmpOpath);
        end
        
    catch ME
        ME
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

function [] = pageTrainingData()

end

function [] = simplePageTrainingData(nX,nY)
    for sel = 1:numel(nX)
        whenPre = whenPoppedNET.predict(nX{sel});
        doesPre = isPoppedNET.predict(nX{sel});
        close all
        plot(nX{sel}','k');
        if ~isempty(doesPre);doesPop = doesPre(1)>.5;else;doesPOP = NaN;end
        title(num2str(doesPop))
        hold on
        plot((whenPre(2,:)>.5) +2,'r--');
        plot(double(nY{sel})-1+1,'b');
        legend({'Predicted','Hand Scored'})
        whenPoppedNET.resetState;
        waitforbuttonpress
        drawnow
    end
end


function [tmpColorSample] = sampleMaskOnly(image,kidx)
    tmpColorSample = [];
    for k = 1:3
        t = double(image(:,:,k));
        tmpColorSample = [tmpColorSample ; t(kidx)];
    end
end

function [tmpColorSample] = sampleRings(image,maskIDX)
    tmpColorSample = [];
    for k = 1:3
        t = double(image(:,:,k));
        featureVec = [];
        for r = 1:numel(maskIDX)
            vec = t(maskIDX{r});
            featureVec = [featureVec ; mean(vec);std(vec,1,1)];
        end
        tmpColorSample = [tmpColorSample ; featureVec];
    end
end


%{
    masterPath = '/mnt/snapper/nate/redCapTrainingData/Emergence_20190519-2019-09-05-15-08-48.4/';
    masterPath = '/mnt/snapper/nate/redCapTrainingData/handTraining/';
    learnRedCappPopper(masterPath)
%}


%{

        % green-ness measurement
        LABsig = {};
        LABsig2 = {};
        parfor maxCHOOSE = 1:100
            gidx = kp(kidx(deltaFrame(maxCHOOSE)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scan and load the data
            mFilePath = nLP{gidx};
            mFileList = {};
            FileExt = {'tif'};
            mFileList = sdig(mFilePath,mFileList,FileExt,1);
            mFileList = mFileList{1};
            tkp = [];
            for e = 1:numel(mFileList)
                tkp(e) = contains(mFileList{e},'raw');
            end
            mFileList = mFileList(logical(tkp));
            NM = [];
            for e = 1:numel(mFileList)
                [pth,nm,ext] = fileparts(mFileList{e});
                fidx = strfind(nm,'_');
                nm = nm(1:(fidx(1)-1));
                NM(e) = str2num(nm);
            end
            [~,sidx] = sort(NM);
            mFileList = mFileList(sidx);
            
            tLABsig = [];
            tLABsig2 = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            for e = 1:numel(mFileList)
                tmpI = double(imread(mFileList{e}))/255;
                tmpI = imresize(tmpI,imageScale2);
                tmpI = bsxfun(@times,tmpI,centerMask);
                tmpLab = rgb2lab(tmpI);
                imshow(tmpLab,[]);
                drawnow
                szLAB = size(tmpLab);
                tLABsig(e,:) = mean(reshape(tmpLab,[prod(szLAB(1:2)) szLAB(3)]),1);
                tLABsig2(e,:) = std(reshape(tmpLab,[prod(szLAB(1:2)) szLAB(3)]),1,1);
                e
            end
            LABsig{maxCHOOSE}(:,:) = tLABsig;
            LABsig2{maxCHOOSE}(:,:) = tLABsig2;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot green-ness
        close all
        toAVG = 50;
        dropThresh = .5;
        greenThresh = 1;
        for e = 1:numel(LABsig)
            
            sig = LABsig{e}(:,2);
            sig2 = LABsig2{e};
            sig = imfilter(sig,fspecial('average',[11 1]),'replicate');
            baseLine = mean(sig(1:toAVG));
            sig = sig - baseLine;
            greenPulse = sig < -dropThresh;
            greenDELTA = min(sig);
            isGreen = greenDELTA < -greenThresh;
            
            
            
            % get index for dataset
            gidx = kp(kidx(deltaFrame(e)));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scan and load the data
            mFilePath = nLP{gidx};
            mFileList = {};
            FileExt = {'tif'};
            mFileList = sdig(mFilePath,mFileList,FileExt,1);
            mFileList = mFileList{1};
            tkp = [];
            for f = 1:numel(mFileList)
                tkp(f) = contains(mFileList{f},'raw');
            end
            mFileList = mFileList(logical(tkp));
            NM = [];
            for f = 1:numel(mFileList)
                [pth,nm,ext] = fileparts(mFileList{f});
                fidx = strfind(nm,'_');
                nm = nm(1:(fidx(1)-1));
                NM(f) = str2num(nm);
            end
            [~,sidx] = sort(NM);
            mFileList = mFileList(sidx);
            
            tmpI = imread(mFileList{end});
            
            
            gridx = find(greenPulse);
            if ~isempty(gridx) & isGreen
                tmpIG = imread(mFileList{gridx(1)});
                CL = 'g';
            else
                tmpIG = tmpI;
                CL = 'k';
            end
            
            
            figure;
            subplot(1,2,1)
            imshow(cat(1,tmpI,tmpIG),[]);
            
            
            subplot(1,2,2)
            plot(sig,'r');
            hold on
            plot(sig2,'b')
            plot(-greenPulse*5,CL)
            %axis([0 size(sig,1) -10 2])
            drawnow
            waitforbuttonpress
            
            close all
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        





        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view results
       

        % scan and show scatter plot
        realValue = [];
        predictValue = [];
        for sel = 1:numel(kp) 
            whenPoppedNET = whenPoppedNET.resetState;
            whenPre = whenPoppedNET.predict(nX{kp(sel)});
            sig = whenPre(2,:) > .5;
            %sig = bwlarge(sig);
            preY = find(sig);
            %{
            
            if whenPre(3,end) > .5
                sigF = whenPre(3,:) > .5;
                sigF = imdilate(sigF,strel('disk',51)) == 1;
                % find movment locations
                pks = find((imdilate(whenPre(2,:),strel('disk',15)) == whenPre(2,:)) & sigF);
                preY = pks(1);
            end
            %}
            %doesPre = isPoppedNET.predict(nX{kp(sel)});
            %doesPre = whenPre(:,end);
            tmpY = find(double(nY{kp(sel)})==2);
            if ~isempty(tmpY)
                realValue(sel) = tmpY(1);
            else
                realValue(sel) = 0;
            end
            if ~isempty(preY)
                predictValue(sel) = preY(1);
            else
                predictValue(sel) = 0;
            end
            
            
            
            if abs(predictValue(sel) - realValue(sel)) > 10
              %  break
            end
            fprintf(['Done with:' num2str(sel) ':' num2str(numel(kp)) '\n']);
        end

        
        
        maxCHOOSE = 1;
        
        kidx = find(realValue ~= 0 & predictValue ~= 0);
        kidx = 1:numel(realValue);
        plot(realValue(kidx),predictValue(kidx),'b.')
        
        [badValues,deltaFrame] = sort(abs(predictValue(kidx) - realValue(kidx)),'descend');
        hold on
        plot(realValue(kidx(deltaFrame(maxCHOOSE))),predictValue(kidx(deltaFrame(maxCHOOSE))),'ro')
        
        
        
        for maxCHOOSE = 1:100

            nLP(kp(kidx(deltaFrame(maxCHOOSE))))
            nY{kp(kidx(deltaFrame(maxCHOOSE)))};
            testPre = whenPoppedNET.predict(nX{kp(kidx(deltaFrame(maxCHOOSE)))});


            sig = testPre(2,:) > .5;
            %sig = bwlarge(sig);

            frP = find(sig > .5);
            frP(1);

            predictValue(kidx(deltaFrame(maxCHOOSE)));
            realValue(kidx(deltaFrame(maxCHOOSE)));
            kp(kidx(deltaFrame(maxCHOOSE)));
            
            testPre(1);
            close all
            plot(nX{kp(kidx(deltaFrame(maxCHOOSE)))}');
            title(num2str(doesPre(1) > .5))
            hold on
            plot(sig +2,'r--');
            hold on
            plot(double(nY{kp(kidx(deltaFrame(maxCHOOSE)))})-1+1)
            whenPoppedNET.resetState;
            drawnow
            waitforbuttonpress
            A = questdlg('Movie');
            if strcmp(A,'Yes')
                mFilePath = nLP{kp(kidx(deltaFrame(maxCHOOSE)))};
                mFileList = {};
                FileExt = {'tif'};
                mFileList = sdig(mFilePath,mFileList,FileExt,1);
                mFileList = mFileList{1};
                tkp = [];
                for e = 1:numel(mFileList)
                    tkp(e) = contains(mFileList{e},'raw');
                end
                mFileList = mFileList(logical(tkp));
                NM = [];
                for e = 1:numel(mFileList)
                    [pth,nm,ext] = fileparts(mFileList{e});
                    fidx = strfind(nm,'_');
                    nm = nm(1:(fidx(1)-1));
                    NM(e) = str2num(nm);
                end
                [~,sidx] = sort(NM);
                mFileList = mFileList(sidx);
                
                handF = find(double(nY{kp(kidx(deltaFrame(maxCHOOSE)))}) == 2);
                if ~isempty(handF)
                    handF = handF(1);
                else
                    handF = NaN;
                end
                
                for e = 1:numel(mFileList)
                    I = imread(mFileList{e});
                    imshow(I,[]);
                    title([num2str(e) '-' num2str(handF)])
                    if e == handF | e == frP(1)
                        waitforbuttonpress
                    end
                    drawnow
                end
               
            end
        end
        
        
        
        
        
        
        maxCHOOSE = 1;
        mFilePath = nLP{kp(kidx(deltaFrame(maxCHOOSE)))};
        mFileList = {};
        FileExt = {'tif'};
        mFileList = sdig(mFilePath,mFileList,FileExt,1);
        mFileList = mFileList{1};
        tkp = [];
        for e = 1:numel(mFileList)
            tkp(e) = contains(mFileList{e},'raw');
        end
        mFileList = mFileList(logical(tkp));
        NM = [];
        for e = 1:numel(mFileList)
            [pth,nm,ext] = fileparts(mFileList{e});
            fidx = strfind(nm,'_');
            nm = nm(1:(fidx(1)-1));
            NM(e) = str2num(nm);
        end
        [~,sidx] = sort(NM);
        mFileList = mFileList(sidx);
        handF = find(double(nY{kp(kidx(deltaFrame(maxCHOOSE)))}) == 2);
        handF = handF(1);
        for e = 1:numel(mFileList)
            I = imread(mFileList{e});
            imshow(I,[]);
            title([num2str(e) '-' num2str(handF)])
            if e == handF
                waitforbuttonpress
            end
            drawnow
        end
        
        
        
        
        
        
        
        sidx = find(predictValue(kidx) == 24);
        predictValue(sidx)
        whenPoppedNET = whenPoppedNET.resetState;
        testPre = whenPoppedNET.predict(nX{kp((sidx(1)))});
        testPre = find(testPre(2,:) > .5)
        LP(kp((sidx(1)))
        plot(nX{kp((sidx(1)))}')
        
        realValue(kp((sidx(1))))
        

%}