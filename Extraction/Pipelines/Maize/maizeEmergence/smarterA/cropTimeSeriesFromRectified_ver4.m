function [finalScore,nX,rawX] = cropTimeSeriesFromRectified_ver4(FileList,boundingBox,timeSeries,tform,cameraShift,emergenceNet,frameNetCourse,frameNetFine,nsz,zU1,zE1,zL1,zERR1,zU2,zE2,zL2,zERR2,kidx)
    % create date: March 18, 2020
    % removed Tscore,LABELS as input
    % FileList      := fileName list of images in the stack
    % boundingBox   := list of boundindBoxs [N by 4]
    % timeSeries    := number sequence for addressing the image list
    resizeAmount = .75;
    
    funcSample = @(X)sampleMaskOnly(X,kidx);
    
    try
        
        
        % removed movVIEW
        %movVIEW = false;
        movSTORE = 1;
        score = zeros(numel(boundingBox),2,numel(timeSeries));
        scoreM = score;
        
        
        cropBoxList = 1:numel(boundingBox);
       
        
        keepMovie = true;
        if keepMovie
            cropBoxList = [4,10,15,19,27,39,57];
            %singleSize = size(tmpI);
            singleSize = [301 301 3];
            movieStore = zeros([singleSize numel(cropBoxList) (numel(timeSeries)-1)]); 
            movieSize = size(movieStore);
        end
        
        
        simMovie = true;
        if simMovie
            simStore1 = zeros([singleSize numel(cropBoxList) (numel(timeSeries)-1)]); 
            simStore2 = zeros([singleSize numel(cropBoxList) (numel(timeSeries)-1)]); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % preallocate matric for pca score data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        PCA_scores = zeros(size(zE1,2)+size(zE2,2),numel(cropBoxList),numel(timeSeries)-1);
        ERROR_scores = zeros(2,numel(cropBoxList),numel(timeSeries)-1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % these look non-used - March 18, 2020
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %spaceFunc{1} = @(X)(X.*bwlarge(X > .0001));
        %spaceFunc{2} = @(X)X;
        %spaceFunc{3} = @(X)X;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each frame in the time series
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        parfor t = 1:(numel(timeSeries)-1)
            fprintf(['*******************************************\n'])
            fprintf(['starting crop on image frame:' num2str(t) ':' num2str(timeSeries(t)) ':' num2str(numel(timeSeries)) '\n']);tic
            fprintf(['*******************************************\n'])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % in the parent call - this is hardwired to empty
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(cameraShift)
                cShift = '';
            else
                cShift = cameraShift{t+1};
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get rectified image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            try
                tmp = getRectifiedImage(FileList{timeSeries(t+1)},tform,[]);
                tmp = tmp * 255;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get rectified image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % these are not used? - march 18, 2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %tmpScore = zeros(numel(boundingBox),2);
            %tmpScoreM = tmpScore;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % preallocate the scores for this frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%]
            tmpPCA_scores = zeros(size(zE1,2)+size(zE2,2),numel(cropBoxList));
            tmpERROR_scores = zeros(2,numel(numel(cropBoxList)));
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % these are not used? - march 18, 2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %L = [];
            %storePath = '';
            
            
            if keepMovie
                tmpStore = zeros(movieSize(1:4));
            end
            
            if simMovie
                tmpSIM1 = zeros(movieSize(1:4));
                tmpSIM2 = zeros(movieSize(1:4));
            end
            
            for cb = 1:numel(cropBoxList)
                fprintf(['starting crop box: ' num2str(cb) ' : ' num2str(numel(cropBoxList)) '\n']);
                
                
                
                
                % crop the cb-th cell from the rectified image
                tmpI = imcrop(tmp,boundingBox{cropBoxList(cb)});
                
                oSZ = size(tmpI);
                if keepMovie
                    try
                        tmpStore(:,:,:,cb) = tmpI;
                    catch
                        fprintf(['ERROR!***************\n']);
                        fprintf([num2str(t) ':' num2str(cb) '-->' size(tmpI) '\n'])
                        fprintf(['ERROR!***************\n']);
                    end
                end
                
                
                % resize image
                tmpI = imresize(tmpI,resizeAmount);
                % sample imge
                tmpSample = funcSample(tmpI);
                
                
                % project image
                cc1 = PCA_REPROJ_T(tmpSample,zE1,zU1);
                simc1 = PCA_BKPROJ_T(cc1,zE1,zU1);
                % error on image
                errorL1 = mean(sum((tmpSample - simc1).*(tmpSample - simc1),1));
                
                
                % project image
                cc2 = PCA_REPROJ_T(tmpSample,zE2,zU2);
                simc2 = PCA_BKPROJ_T(cc2,zE2,zU2);
                % error on image
                errorL2 = mean(sum((tmpSample - simc2).*(tmpSample - simc2),1));
                
                
                 if simMovie
                    try
                        tmpS = putDataIntoImageAtP(simc1,kidx,size(tmpI));
                        tmpS = imresize(tmpS,oSZ(1:2));
                        tmpSIM1(:,:,:,cb) = tmpS;
                        
                        tmpS = putDataIntoImageAtP(simc2,kidx,size(tmpI));
                        tmpS = imresize(tmpS,oSZ(1:2));
                        tmpSIM2(:,:,:,cb) = tmpS;
                    
                    catch
                        fprintf(['ERROR!***************\n']);
                        fprintf([num2str(t) ':' num2str(cb) '-->' size(tmpI) '\n'])
                        fprintf(['ERROR!***************\n']);
                    end
                end
                
                %{
                tmpM = zeros(size(tmpI,1),size(tmpI,2));
                tmpM(kidx) = 1;
                tmpV = reshape(bindVec(tmpI(:)),size(tmpI));
                out = flattenMaskOverlay(tmpV,logical(tmpM),.7,'b');
                imshow(out,[]);
                drawnow
                %}
                tmpPCA_scores(:,cb) = [cc1;cc2];
                tmpERROR_scores(:,cb) = [errorL1;errorL2];
               
                fprintf(['image cropped of size:[' num2str(size(tmpI)) ']\n']);
                %[pth,nm,ext] = fileparts(FileList{timeSeries(t)});
                
                %{
                if ~isempty(CELLMASK)
                    
                    %VV = zeros(3*numel(vidx),1);
                   
                    
                    tmpI = bsxfun(@times,255*single(imresize(tmpI,[nsz,nsz])),CELLMASK);
                    
                    %tmpI = bsxfun(@times,255*single(imresize(tmpI,[nsz,nsz])),CELLMASK);
                    %%%%%%%%%%%%%%%%
                    % the line below might be causing errors
                    % I found the errors a while ago 
                    % This line is being put back in on Nov 21, 2019
                    tmpI = imcrop(tmpI,BOX);
                    %%%%%%%%%%%%%%%%
                    %{
                    if movSTORE
                        temp2(:,:,:,t) = tmpI;
                    end
                    %}
                    
                    
                    
                    VV = [];
                    for k = 1:size(tmpI,3)
                        tmpSS = tmpI(:,:,k);
                        VV = [VV ;tmpSS(kidx)];
                    end
                    
                    
                    %tmpI = rgb2hsv_fast(tmpI/255);
                    %tmpI(:,:,3) = imhistmatch(tmpI(:,:,3),nI);
                    %tmpI = 255*hsv2rgb(tmpI);
                    
                    %if movSTORE
                    % remove the store -  i think this was for saving - Nov
                    % 21 2019
                    % 
                    %temp(:,:,:,t) = tmpI;
                    %end
                    
                    %tmpI = tmpI - imcrop(single(imresize(baseline{cb},[nsz,nsz])),subBOX);
                end
                %}
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % remove too?
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %{
                if movVIEW
                    %MM = bsxfun(@times,255*single(imresize(tmpI,[nsz,nsz])),CELLMASK);
                    imshow(tmpI/255,[]);
                    STACKVIEW(:,:,:,vc) = tmpI/255;
                    vc = vc + 1;
                    title(num2str(t))
                    drawnow
                end
                %}
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %{
                % BLOCK for movVIEW - removed march 18, 2020
                if ~movVIEW
                    

                   

                    if isempty(emergenceNet)
                        %tmpMODEL = PCA_REPROJ_T(tmpI(:),zE,zU);
                        %tmpMODEL = PCA_BKPROJ_T(tmpMODEL,zE,zU);
                        
                        
                        %sz = size(tmpI);
                        %tmpIM = reshape(tmpMODEL,sz);
                        
                        %{
                        % make disk mask
                        [idx,model,POS,out] = getCellModel2(tmpIM,fullLABEL{1},fullLABEL{2},1,spaceFunc,false);
                        % copy over raw data into model for fusion
                        repidx = find(bwlarge(idx==1));
                        for k = 1:3
                            tmpp = tmpI(:,:,k);
                            tmpTar = tmpIM(:,:,k);
                            tmpTar(repidx) = tmpp(repidx);
                            tmpIM(:,:,k) = tmpTar;
                        end
                        [mea,POSRED] = splitEnergyLevels(tmpIM,idx,HI,bwlarge(idx==1));
                        
                        
                        rPOS = POS;
                        tL = [];
                        for k = 1:3
                            tL(k) = norm(rPOS(:,1)-rPOS(:,2));
                            %plot(rPOS(2,1:2),rPOS(1,1:2),'k');
                            rPOS = circshift(rPOS,[0 1]);
                        end
                        tL = [tL' ; norm(POS(:,1) - POSRED')];
                        L(:,cb) = tL;
                        %}
                        %%{
                        
                        
                        %{
                            CL = {'r' 'c' 'y'};
                            RGB = label2rgb(idx);
                            if disp
                                figure(h1);
                                %imshow(cat(2,tmpI/255,out,double(RGB)/255),[]);
                                imshow(tmpIM/255,[]);
                                hold on
                                for k = 1:size(POS,2)
                                    plot(POS(2,k),POS(1,k),[CL{k} 'o']);
                                    hold on
                                end
                                plot(POSRED(1,2),POSRED(1,1),'r*');

                                rPOS = POS;
                                for k = 1:3
                                    plot(rPOS(2,1:2),rPOS(1,1:2),'k');
                                    rPOS = circshift(rPOS,[0 1]);
                                end
                                plot([POS(2,1) POSRED(2)],[POS(1,1) POSRED(1)],'w');
                                plot(size(RGB,2)/2,size(RGB,1)/2,'k.');
                                hold off
                                title(num2str(e))
                                sL = [];
                                %kL = squeeze(LENM(:,cb,:));
                                
                                for s = 1:size(kL,1)
                                    sL(s,:) = imfilter(kL(s,:),fspecial('average',[1 10]),'replicate');
                                end
                                
                                TT = 10;
    
                                initL = mean(sL(:,1:(min(t-1,TT))),2);
                                STR = bsxfun(@minus,sL,initL);
                                STR = bsxfun(@minus,STR,initL.^-1);
                                STR(end+1,:) = mean(abs(STR),1);
                                figure(h2);
                                plot(STR');
                                hold on
                                %plot(TYPEF(e,1:t)*5,'k')
                                hold off
                                drawnow
                                e
                            end
                        %%}
                        
                        
                        %{
                        tmpMODEL = reshape(tmpMODEL,size(tmpI));
                        tmpScore(cb,:) = predict(emergenceNet,tmpI);
                        tmpScoreM(cb,:) = predict(emergenceNet,tmpMODEL);
                        %}
                        
                        
                        
                        
                        
                        if movSTORE
                            temp(:,:,:,t) = tmpIM;
                        end
                        %}
                    end

                    if (~isempty(zU) && ~isempty(zE))
                        %tmpPCA_scores(:,cb) = PCA_REPROJ_T(tmpI(:),zE,zU);
                        tmpPCA_scores(:,cb) = PCA_REPROJ_T(VV,zE,zU);
                    end


                    %{
                    %%%% removing the spool to disk option as of Nov 21, 2019
                    % store the tmpI image under the key
                    if ~isempty(storePath)
                        fidx = strfind(pth,filesep);
                        key = [storePath filesep lower((pth((fidx(end)+1):end))) '-' LABELS{cb} '--f' num2str(t) '.tif'];
                        imwrite(tmpI,key);
                    end
                    %}


                    fprintf(['ending crop box: ' num2str(cb) ' : ' num2str(numel(boundingBox)) '\n']);
                end
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
            end
            if keepMovie
                movieStore(:,:,:,:,t) = tmpStore;
            end
            
            if simMovie
                simStore1(:,:,:,:,t) = tmpSIM1;
                simStore2(:,:,:,:,t) = tmpSIM2;
            end
            
            %{
            %%%% This looks like it is needed and no longer should be 
            %%%% under the flag of movVIEWas of Nov 21, 2019
            %}
            %if ~movVIEW
                %LENM(:,:,t) = L;
                %score(:,:,t) = tmpScore;
                %scoreM(:,:,t) = tmpScoreM;
                
                
               
                    PCA_scores(:,:,t) = tmpPCA_scores;
                    ERROR_scores(:,:,t) = tmpERROR_scores;
                
            %end

            fprintf(['**************************************************************************************\n'])
            fprintf(['ending crop on image frame:' num2str(toc) ':'  num2str(t) ':' num2str(timeSeries(t)) ':' num2str(numel(timeSeries)) '\n']);
            fprintf(['**************************************************************************************\n'])
        end
        
        
        %
        %cropBoxList = [4,10,15,19,27*,39,57];
        toView = 5;
        for e = 1:(finePrediction(toView)+20)%size(movieStore,5)
            
            imshow([movieStore(:,:,:,toView,e) simStore1(:,:,:,toView,e) simStore2(:,:,:,toView,e)]/255,[]);
            %{
            title([num2str(e) '--' num2str(cropBoxList(toView)) '--' ...
                num2str(finePrediction(cropBoxList(toView)))]);
  
            if e == finePrediction(cropBoxList(toView))
                waitforbuttonpress
            end
            %}
            
            title([num2str(e) '--' num2str(cropBoxList(toView)) '--' ...
                num2str(finePrediction(toView))]);
  
            if e == finePrediction(toView)
                waitforbuttonpress
                
            end
            
            drawnow
        end
        %
        
        
        W = 15;
        numCom = 7;
        toUse = [[1:7],[10:17]];
        for cell = 1:size(PCA_scores,2)
            fprintf(['Grading cell:' num2str(cell) ':' num2str(size(PCA_scores,2)) '\n']);
            sig = squeeze(PCA_scores(:,cell,:));
            sig_error = squeeze(ERROR_scores(1:2,cell,:));
            sig = [sig;sig_error];
            sig = bsxfun(@minus,sig,sig(:,1));
            % eigenvalue normalize - with error
            sig = bsxfun(@times,sig,[zL1;zL2;zERR1;zERR2].^-.5);
            sig = sig(toUse,:);
            tmp = double(emergenceNet.predict(sig));
            emergenceNet = emergenceNet.resetState;
            m1(cell) = tmp(2);
            tmp = frameNetCourse.predict(sig);
            frameNetCourse = frameNetCourse.resetState;
            m2(cell,:) = tmp(2,:);
            bsig = m2(cell,:) > .5;
            bidx = find(bsig);
            
            coursePrediction(cell) = 0;
            if ~isempty(bidx);coursePrediction(cell) = bidx(1);end
            
            finePrediction(cell) = coursePrediction(cell);
            
            if coursePrediction(cell) ~= 0
                str = bidx(1);
                str = str - W;
                if str <= 0;str = 16;end
                if str >= (size(sig,2) - 2*W);str = (size(sig,2) - 2*W);end
                widx = str:(str+2*W);
                wsig = sig(:,widx);
                tmp = frameNetFine.predict(wsig);
                frameNetFine = frameNetFine.resetState();
                m3(cell,:) = tmp(2,:);
                bsig = m3(cell,:) > .5;
                % offset back by W
                bidx = find(bsig) - (W + 1);
                if ~isempty(bidx);finePrediction(cell) = coursePrediction(cell) + bidx(1);end
            end
        end
        
        LONG_term = 1;
        NUMPC = 10;
        PCA_scores = PCA_scores(1:NUMPC,:,:);
        
        
        
        
        
        % for quantum
        %LENM = permute(LENM,[1 3 2]);
        PCA_scores = permute(PCA_scores,[1 3 2]);
        popThreshold = .6;
        [finalScore,rawX,nX] = gradeRedCaps(PCA_scores,emergenceNet,frameNetCourse,popThreshold,.5);
        %{
        pY = BESTEmergence.predict(nX);
        for e = 1:numel(pY)
            
            [stepSig,finalScore(e)] = extractFrameFromProb(pY{e}(2,:),.4,5,21);
            %{
            [~,midx] = max(,[],1);
            fidx = find(midx==2);
            if ~isempty(fidx)
                finalScore(e) = fidx(1);
            else
                finalScore(e) = 0;
            end
            %}
        end
        %}
        %{
        %scoreN = .5*(score + scoreM);
        fprintf(['Starting image grading.\n']);tic
        if (~isempty(zU) && ~isempty(zE))
            PCA_scores = permute(PCA_scores,[1 3 2]);
            for tr = 1:size(PCA_scores,3)
                PCA_scores(:,:,tr) = myFrenet(PCA_scores(:,:,tr),[7],@(X)funcG(X));
            end
            PCA_scores = abs(PCA_scores);
            % try uniform z normalize
            psz = size(PCA_scores);
            PCA_scores = reshape(PCA_scores,[psz(1) prod(psz(2:3))]);
            PCA_scores = bsxfun(@minus,PCA_scores,Kmu);
            PCA_scores = bsxfun(@times,PCA_scores,Kstd.^-1);
            PCA_scores = reshape(PCA_scores,psz);
            sTHRESH = .5;
            dTHRESH = .5;
            dScore = [];
            sScore = [];
            finalScore = [];
            handScore = [];
            sz = size(PCA_scores);
            zPCA_scores = reshape(PCA_scores,[sz(1) prod(sz(2:3))]);
            zPCA_scores = bsxfun(@minus,zPCA_scores,Zmean);
            zPCA_scores = bsxfun(@times,zPCA_scores,Zstd.^-1);
            zPCA_scores = reshape(zPCA_scores,sz);
            
            
            %zPCA_scores = permute(zPCA_scores,[1 3 2]);
            szZ = size(zPCA_scores);
            zPCA_scores = reshape(zPCA_scores,[szZ(1:2) 1 szZ(3)]);
            wPCA_scores = reshape(PCA_scores,[szZ(1:2) 1 szZ(3)]);

            zPCA_scores = cat(3,zPCA_scores,wPCA_scores,zeros(size(zPCA_scores)));
            
            hSZ = (WINDOW_SZ-1)/2;
            for tr = 1:size(zPCA_scores,4)
                fprintf(['Starting image grading.:' num2str(tr) ':' num2str(size(zPCA_scores,4)) '\n']);
                Mdata = [];
                if ~LONG_term
                    for l = 1:3
                        tmpData = im2col(zPCA_scores(:,:,l,tr),[size(zPCA_scores,1) WINDOW_SZ]);
                        tmpData = reshape(tmpData,[size(zPCA_scores,1) WINDOW_SZ size(tmpData,2)]);
                        tmpSz = size(tmpData);
                        tmpData = reshape(tmpData,[tmpSz(1:2) 1 tmpSz(3)]);
                        Mdata = cat(3,Mdata,tmpData);
                    end
                    
                    wholeP = predict(BESTemergence,Mdata);
                    wholeP = [zeros(hSZ,size(wholeP,2));wholeP;zeros(hSZ,size(wholeP,2))];
                end
                
                if LONG_term
                    %sig = [zPCA_scores(:,:,1,tr);zPCA_scores(:,:,2,tr)];
                    sig = [zPCA_scores(:,:,2,tr)];
                    %sig = bsxfun(@minus,sig,mean(sig(:,1:20),2));
                    %sig = imfilter(sig,fspecial('average',[1 7]),'replicate');
                    [updatedNet,class,wholeP] = classifyAndUpdateState(BESTemergence,{sig});
                    wholeP = wholeP{1}';
                    %wholeP = double(wholeP{1})-1;
                    %wholeP = [wholeP' wholeP'];
                end
                
               
                
                
                dynamicIDX = find(wholeP(:,2) > dTHRESH);
                dProb = wholeP(:,2);
                sig = wholeP(:,2) > dTHRESH;
                %sig = imclose(sig,ones(6,1));
                sig = bwareaopen(sig,12);
                R = regionprops(sig,'PixelIdxList');
                dynamicIDX = zeros(numel(timeSeries),1);
                if ~isempty(R)
                    dynamicIDX(R(1).PixelIdxList) = 1;
                    dynamicIDX = dynamicIDX.*dProb;
                    dynamicIDX = dynamicIDX / sum(dynamicIDX);
                    dynamicScore = dynamicIDX'*timeSeries';
                    dynamicScore = R(1).PixelIdxList(1);
                else
                    dynamicScore = 0;
                end
                
                
                
                
                spaceIDX = squeeze(score(tr,2,:));
                spaceIDX = circshift(spaceIDX,[1 0]);
                sProb = spaceIDX;
                doesGerminate = sum(spaceIDX > sTHRESH) > 7;
                doesGerminate = sum(sig) > 7;
                spaceScore = find(spaceIDX > sTHRESH);
                if ~isempty(spaceScore)
                    spaceScore = spaceScore(1);
                else
                    spaceScore = 0;
                end
                
                
                
                
                spaceIDXM = squeeze(scoreM(tr,2,:));
                spaceIDXM = circshift(spaceIDXM,[1 0]);
                sProb = spaceIDXM;
                %doesGerminate = sum(spaceIDXM > sTHRESH) > 25;
                spaceScoreM = find(spaceIDXM > sTHRESH);
                if ~isempty(spaceScoreM)
                    spaceScoreM = spaceScoreM(1);
                else
                    spaceScoreM = 0;
                end
                
                
                
                dScore(tr) = dynamicScore*doesGerminate;
                sScore(tr) = spaceScore*doesGerminate;
                sScoreM(tr) = spaceScoreM*doesGerminate;
                
                
                
                if doesGerminate && dScore(tr) == 0
                    dProb = imfilter(dProb,fspecial('average',5),'replicate');
                    [~,dScore(tr)] = max(dProb);
                end
                
                
                if dScore(tr) ~= 0 && sScore(tr) ~= 0
                    
                    %%%%
                    toAverage = dScore(tr);
                    if abs(dScore(tr) - sScore(tr)) < 10
                        toAverage = [toAverage sScore(tr)];
                    end
                    
                    if abs(dScore(tr) - sScoreM(tr)) < 10
                        toAverage = [toAverage sScoreM(tr)];
                    end
                    
                    finalScore(tr) = nanmean(toAverage);
                    
                    %%%%%%
                    
                elseif dScore(tr) ~= 0 && sScore(tr) == 0
                    % only dynamic score
                    finalScore(tr) = dScore(tr);
                elseif dScore(tr) == 0 && (sScore(tr) ~= 0 || sScoreM(tr) ~= 0)
                    % dynamic has failed but others have noe
                    toAverage = [];
                    %%%%
                    if sScore(tr) ~= 0 
                        toAverage = [toAverage sScore(tr)];
                    end
                    if sScoreM(tr) ~= 0
                        toAverage = [toAverage sScore(tr)];
                    end
                    
                    finalScore(tr) = nanmean(toAverage);
                    
                    %%%%%%
                    
                elseif dScore(tr) == 0 && sScore(tr) == 0
                    finalScore(tr) = 0;
                end
                
                
                if ~isempty(Tscore)
                    handScore(tr) = Tscore(tr);
                    handScoreProb = zeros(numel(timeSeries),1);
                    if handScore(tr) ~= 0
                        handScoreProb(handScore(tr):end) = 1;
                    end
                end
                
                
                sScoreStep = zeros(numel(timeSeries),1);
                dScoreStep = zeros(numel(timeSeries),1);
                sScoreStepM = zeros(numel(timeSeries),1);
                fScoreStep = zeros(numel(timeSeries),1);
                if sScore(tr) ~= 0
                    sScoreStep(sScore(tr):end) = 1;
                end
                if dScore(tr) ~= 0 
                    dScoreStep(dScore(tr):end) = 1;
                end
                if sScoreM(tr) ~= 0 
                    sScoreStepM(sScoreM(tr):end) = 1;
                end
                if finalScore(tr) ~= 0 
                    fScoreStep(finalScore(tr):end) = 1;
                end
                
               
                
                %{
                plot(1:numel(timeSeries),.8*sScoreStep,'r');
                hold on
                plot(1:numel(timeSeries),.4*dScoreStep,'b');
                if ~isempty(Tscore)
                    plot(1:numel(timeSeries),handScoreProb,'g');
                end
                plot(1:numel(timeSeries),.6*sScoreStepM,'c');
                plot(1:numel(timeSeries),1.2*fScoreStep,'k');
                plot(1:numel(timeSeries),dProb,'c');
                plot(1:numel(timeSeries),sProb,'m');
                plot(1:numel(timeSeries),.5*ones(numel(timeSeries),1),'y')
                axis([0 numel(timeSeries) 0 1.2])
                drawnow
                hold off
                %}
                kidx = handScore ~= 0 & finalScore ~= 0;
                DELTA = mean(abs(handScore(kidx) - finalScore(kidx)));
                %{
                if ~isempty(handScore(kidx))
                    title([num2str(tr) '-' num2str(finalScore(tr)) '-' num2str(handScore(tr)) '----' num2str(corr(finalScore(kidx)',handScore(kidx)')) '------' num2str(DELTA)]);
                else
                    title([num2str(tr) '-' num2str(finalScore(tr)) '-' num2str(handScore(tr)) '----' num2str(0) '------' num2str(DELTA)]);
                end
                drawnow
                %}
                %waitforbuttonpress
                fprintf(['Ending image grading.:' num2str(tr) ':' num2str(size(zPCA_scores,4)) '\n']);
            end
        end
        %}
        fprintf(['Ending image grading: ' num2str(toc) '\n']);
        %{
        for t = 1:numel(tSTACK)
            STACK(:,:,:,:,t) = tSTACK{t};
            tSTACK{t} = [];
        end
        STACK = permute(STACK,[1 2 3 5 4]);
        %}
    catch ME
        getReport(ME)
        here = 1;
    end
        
end



function [tmpColorSample] = sampleMaskOnly(image,kidx)
    tmpColorSample = [];
    for k = 1:3
        t = double(image(:,:,k));
        tmpColorSample = [tmpColorSample ; t(kidx)];
    end
end