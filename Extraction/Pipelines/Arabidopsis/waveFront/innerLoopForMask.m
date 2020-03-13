function [] = innerLoopForMask(STACK,plantMask,baseName,oPath,nm,maskFileList,configData)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INPUTS:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STACK         := the raw-data
        % pM            := the plantMask to use for the processing
        % baseName      := the base file name for the data
        % oPath         := the output location
        % nm            := name of the input file
        % maskFileList  := list of mask files to use for focus
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        baseTHIS = [oPath baseName filesep];
        mkdir(baseTHIS);
        thisOUTPUT = [baseTHIS 'WAVE_ALGO' filesep];
        mkdir(thisOUTPUT);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the configuration data from the struct 
        CL = configData.CL;
        toRenderVideo = configData.toRenderVideo;
        waveThresholdPercent = configData.waveThresholdPercent;
        t_stack_SKIP_frame_vec = configData.t_stack_SKIP_frame_vec;
        PER_contour = configData.PER_contour;
        globalSKIP = configData.globalSKIP;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        analysisMag = 10;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OLD CODE for direct configure parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % colors for displays
        CL = {'r' 'g' 'b'};
        % flag to render video to disk
        toRenderVideo = 1; 
        % which are less than PER of the max
        waveThresholdPercent = [.2 .35 .5];
        % skips for the frame over the t-stacks
        t_stack_SKIP_frame_vec = 1:10:31;
        % precent contour levels to process the data
        PER_contour = [.2 .35 .5];
        % the global skip value for the analysis
        globalSKIP = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract SIG from STACK
        % the signal that is processed is:
        % 1) subtract the initial signal from the stack
        % 2) divide by the max to have each signal eleOf [0,1]
        % 3) if the max == 0, then set to 1
        % 4) filter the signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % subtract each signal to start at zero
        SIG = bsxfun(@minus,STACK,STACK(:,:,1));
        % get the max value
        maxValue = max(SIG,[],3);
        % if max value == 0, set to 1
        maxValue(maxValue==0) = 1;
        % and have each signal go to one as the max
        SIG = bsxfun(@times,SIG,maxValue.^-1);
        sz = size(SIG);
        SIG = reshape(SIG,[prod(sz(1:2)) sz(3)]);
        % filter SIG
        SIG = imfilter(SIG,ones(1,configData.time_filt)/configData.time_filt,'replicate');
        % reshape SIG
        SIG = reshape(SIG,sz);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % renormalize
        % subtract filtered value
        SIG = bsxfun(@minus,SIG,SIG(:,:,1));
        % get the max value
        maxValue = max(SIG,[],3);
        % if max value == 0, set to 1
        maxValue(maxValue==0) = 1;
        % and have each signal go to one as the max
        SIG = bsxfun(@times,SIG,maxValue.^-1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % analysis of each pixel independently
        % find the frame for the transistion into, outof, peak 
        % find the velocity at into and outof
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate the gradient for physics wave velocity
        [d1,d2,d3] = gradient(SIG);
        if ~any(configData.gradient_filt == 0)
            filtBank = fspecial('gaussian',[configData.gradient_filt(1:2)],configData.gradient_filt(3));
            for tm = 1:size(SIG,3)
                d1(:,:,tm) = imfilter(d1(:,:,tm),filtBank,'replicate');
                d2(:,:,tm) = imfilter(d2(:,:,tm),filtBank,'replicate');
            end
        end
        % calculate the physics wave velocity
        velW = d3.*(d1.^2 + d2.^2).^-.5;
        % reshape the velocity
        vSIG = reshape(velW,[prod(sz(1:2)) sz(3)]);
        % reshape the sig
        SIG = reshape(SIG,[prod(sz(1:2)) sz(3)]);
        X = 1:size(SIG,2);
        Xi = linspace(1,size(SIG,2),analysisMag*size(SIG,2));
       
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % model extension(s)
        MODEL_EXT = 50;
        % model threshold
        MOD_THRESH = .2;
        % quick model
        dN = (d1.^2 + d2.^2).^-.5;
        n1 = d1.*dN;
        n2 = d2.*dN;
        Vn1 = reshape(n1,[prod(sz(1:2)) sz(3)]);
        Vn2 = reshape(n2,[prod(sz(1:2)) sz(3)]);
        tmp_n1 = Vn1(p,:);
        tmp_n2 = Vn2(p,:);
        % quick model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % quick model
        frameLook = 100;
        tmpSIG = reshape(SIG,sz);
        
        
        for frameLook = 1:100
            data = double(plantMask.*tmpSIG(:,:,frameLook));
            dataNext = double(plantMask.*tmpSIG(:,:,frameLook+1));
            imshow(data,[0 1]);

            % get the wave front
            dC = contourc(data,[.8 .8]);
            imshow(data,[0 1]);
            hold on
            % init variables for the contours
            cnt = 1;CON = {};
            % while there are still contours
            while ~isempty(dC)
                % get the contour pointer
                ptr = dC(:,1);
                % get the contour
                CON{cnt} = dC(:,2:1+ptr(2));
                % pop/remove the contour from the queue
                dC(:,1:1+ptr(2)) = [];
                % if render is true - render the contour
                plot(CON{cnt}(1,:),CON{cnt}(2,:),'r');hold on


                nii1 = ba_interp2(double(n1(:,:,frameLook)),CON{cnt}(1,:),CON{cnt}(2,:));
                nii2 = ba_interp2(double(n2(:,:,frameLook)),CON{cnt}(1,:),CON{cnt}(2,:));
                speed = ba_interp2(double(velW(:,:,frameLook)),CON{cnt}(1,:),CON{cnt}(2,:));

                
                
                if any(speed(:) > 10)
                    frameLook;
                end
                
                
                plot(CON{cnt}(1,:)-nii1.*speed,CON{cnt}(2,:)-nii2.*speed,'m')
                


                %quiver(CON{cnt}(1,:),CON{cnt}(2,:),nii1.*speed,nii2.*speed,'Color','r');
                % increment the contour count
                cnt = cnt + 1;
            end


            % get the wave front
            dC = contourc(dataNext,[.8 .8]);
            % init variables for the contours
            cnt = 1;CON = {};
            % while there are still contours
            while ~isempty(dC)
                % get the contour pointer
                ptr = dC(:,1);
                % get the contour
                CON{cnt} = dC(:,2:1+ptr(2));
                % pop/remove the contour from the queue
                dC(:,1:1+ptr(2)) = [];
                % if render is true - render the contour
                plot(CON{cnt}(1,:),CON{cnt}(2,:),'b');hold on


                %nii1 = ba_interp2(double(n1(:,:,frameLook)),CON{cnt}(1,:),CON{cnt}(2,:));
                %nii2 = ba_interp2(double(n2(:,:,frameLook)),CON{cnt}(1,:),CON{cnt}(2,:));
                %speed = ba_interp2(double(velW(:,:,frameLook)),CON{cnt}(1,:),CON{cnt}(2,:));


                %quiver(CON{cnt}(1,:),CON{cnt}(2,:),nii1.*speed,nii2.*speed);
                % increment the contour count
                cnt = cnt + 1;
            end
            hold off
            pause(.1);
            if frameLook == 20
                %waitforbuttonpress
            end
        end
        
        
        % quick model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over each pixel
        parfor e = 1:size(SIG,1)
            % interpolate the raw signal
            tmpNormalizedSignal = interp1(X,SIG(e,:),Xi,'spline');
            % interpolate the velocity signal
            tmpVelocitySignal = interp1(X,vSIG(e,:),Xi,'linear');
            
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % quick model
            tmpNormalizedSignal_1 = interp1(X,Vn1(e,:),Xi,'spline');
            tmpNormalizedSignal_2 = interp1(X,Vn2(e,:),Xi,'spline');
            % quick model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            
            % create a bank of variables
            TMPpeakFrame = zeros(1,numel(waveThresholdPercent));
            TMPinFrame = zeros(1,numel(waveThresholdPercent));
            TMPvelFrontFrame = zeros(1,numel(waveThresholdPercent));
            TMPoutFrame = zeros(1,numel(waveThresholdPercent));
            TMPvelBackFrame = zeros(1,numel(waveThresholdPercent));
            
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % quick model
            try
                modelBorder = [];
                modelRegion = tmpNormalizedSignal > MOD_THRESH;
                modelRegion = bwlarge(modelRegion);
                midx = find(modelRegion);
                [~,maxIDX] = max(tmpNormalizedSignal);
                
                riseTM = maxIDX - midx(1);
                recoverTM = midx(1) - maxIDX;
                riseSpeed = tmpVelocitySignal(midx(1));
                recoverSpeed = tmpVelocitySignal(midx(end));
                
                modelBorder(1) = max(midx(1)-MODEL_EXT,1);
                modelBorder(2) = maxIDX;
                modelBorder(3) = min(midx(end)+MODEL_EXT,numel(modelRegion));
                modelWaveRise = tmpNormalizedSignal(modelBorder(1):modelBorder(2));
                modelWaveRecover = tmpNormalizedSignal((modelBorder(2)+1):modelBorder(3));
                ImodelWaveRise = interp1(1:numel(modelWaveRise),modelWaveRise,linspace(1,numel(modelWaveRise),200));
                ImodelWaveRecover = interp1(1:numel(modelWaveRecover),modelWaveRecover,linspace(1,numel(modelWaveRecover),200));

                phase1(e,:) = ImodelWaveRise;
                phase2(e,:) = ImodelWaveRecover;
                waveFeatures(e,:) = [riseTM recoverTM riseSpeed recoverSpeed];
                modelDomain(e,:) = modelBorder;
                modelUse(e) = true;
            catch
                phase1(e,:) = zeros(1,200);
                phase2(e,:) = zeros(1,200);
                waveFeatures(e,:) = zeros(1,4);
                modelDomain(e,:) = zeros(1,3);
                modelUse(e) = false;
            end
        end
            
        %}
            
        %{
        
            % get the indexing and remove the "bad" points
            pidx = find(plantMask);
            fidx = find(modelUse);
            pidx = intersect(pidx,fidx);
            [mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(double([phase1(pidx,:) phase2(pidx,:)]),3);
            plot(mC(:,1),mC(:,2),'.')
            kp = abs(mC(:,1)) < 30 & abs(mC(:,2)) < 5;
            kp = find(kp);
            
            
            
            % model together
            [mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(double([phase1(pidx(kp),:) phase2(pidx(kp),:)]),3);
            whenPeak = modelDomain(fidx(kp),2);
            plot(mC(:,1),mC(:,2),'.')
            fidx = find(modelUse);
            toView = 1000;
            plot([phase1(pidx(1:toView),:) phase2(pidx(1:toView),:)]')
            for e = 1:size(mS,1)
                plot(mS(e,:),'r');
                hold on
                plot([phase1(pidx(kp(e)),:) phase2(pidx(kp(e)),:)],'k');
                hold off 
                pause(.1)
            end
            % model seperate
            [mS1 mC1 mU1 mE1 mL1 mERR1 mLAM1] = PCA_FIT_FULL(double([phase1(pidx(kp),:)]),3);
            [mS2 mC2 mU2 mE2 mL2 mERR2 mLAM2] = PCA_FIT_FULL(double([phase2(pidx(kp),:)]),3);
            whenPeak = modelDomain(fidx(kp),2);
            for e = 1:size(mS,1)
                plot([mS1(e,:) mS2(e,:)],'r');
                hold on
                plot([phase1(pidx(kp(e)),:) phase2(pidx(kp(e)),:)],'k');
                hold off 
                pause(.1)
            end
            
            
            % look at temporal distribution of model parameters
            close all
            toMake = 1;
            modelImage = zeros(size(plantMask));
            modelImage(pidx(kp)) = mC1(:,toMake);
            range = [min(mC1(:,toMake)) max(mC1(:,toMake))];
            modelImage = reshape(modelImage,size(plantMask));
            imshow(modelImage,range);
            waitforbuttonpress
            close all
            toMake = 2;
            modelImage = zeros(size(plantMask));
            modelImage(pidx(kp)) = mC1(:,toMake);
            range = [min(mC1(:,toMake)) max(mC1(:,toMake))];
            modelImage = reshape(modelImage,size(plantMask));
            imshow(modelImage,range);
            waitforbuttonpress
            close all
            toMake = 3;
            modelImage = zeros(size(plantMask));
            modelImage(pidx(kp)) = mC1(:,toMake);
            range = [min(mC1(:,toMake)) max(mC1(:,toMake))];
            modelImage = reshape(modelImage,size(plantMask));
            imshow(modelImage,range);
            
            
            % make simulated movie
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over each pixel
        mM = zeros(size(SIG));
        pidx2 = pidx(kp);
        parfor e = 1:numel(pidx2)
            sig = zeros(size(Xi));
            sigX = [mS1((e),:) mS2((e),:)];
            sigX = interp1(1:numel(sigX),sigX,linspace(1,numel(sigX),modelDomain(pidx2(e),3)-modelDomain(pidx2(e),1)+1));
            sig(modelDomain(pidx2(e),1):modelDomain(pidx2(e),3)) = sigX;
            sig = interp1(1:numel(sig),sig,linspace(1,numel(sig),size(SIG,2)));
            BANK(e,:) = sig;
            LOC(e) = pidx2(e);
            %mM(pidx(e),:) = sig;
            e
            numel(pidx2)
            
        end
        mM(LOC,:) = BANK;
        for e = 1:size(mM,3)
            imshow(mM(:,:,e),[0 1]);
            drawnow
        end
            
            % link between rise and recover
            r1 = 1;
            r2 = 1;
            close all
            plot(mC1(:,r1),mC2(:,r2),'.')
            
            % parameter sweep
            STEPS = 10;
            for toSweep = 1:size(mC1,2)
                [sweepD] = sweepPCA(mC1,mE1,mU1,diag(mLAM1).^.5,toSweep,STEPS);
                for e = 1:size(sweepD,2)
                    plot(squeeze(sweepD(1,e,:)))
                    hold on
                end
                waitforbuttonpress
                close all
            end
            
            
               
            % parameter sweep
            STEPS = 10;
            for toSweep = 1:size(mC2,2)
                [sweepD] = sweepPCA(mC2,mE2,mU2,diag(mLAM2).^.5,toSweep,STEPS);
                for e = 1:size(sweepD,2)
                    plot(squeeze(sweepD(1,e,:)))
                    hold on
                end
                waitforbuttonpress
                close all
            end
            
            
            % look at plots of wave features vs pca scores
            close all
            plot(whenPeak,mC1(:,1),'.')
            waitforbuttonpress
            close all
            plot(waveFeatures(pidx(kp),1),mC1(:,1),'.')
            waitforbuttonpress
            close all
            plot(waveFeatures(pidx(kp),2),mC1(:,1),'.')
             waitforbuttonpress
            close all
            plot(waveFeatures(pidx(kp),3),mC1(:,1),'.')
            waitforbuttonpress
            close all
            plot(waveFeatures(pidx(kp),4),mC1(:,2),'.')
            
            
            % quick model
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
            %}
        
    
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for each wave threshold - perform the analysis
            for per = 1:numel(waveThresholdPercent)
                % binary vector for wave transistions
                waveTransistionINTO = tmpNormalizedSignal(1:end-1) < waveThresholdPercent(per) & tmpNormalizedSignal(2:end) > waveThresholdPercent(per);
                waveTransistionOUTOF = tmpNormalizedSignal(1:end-1) > waveThresholdPercent(per) & tmpNormalizedSignal(2:end) < waveThresholdPercent(per);

                % get the index of the wave transistions and the wave peak
                idxINTO = find(waveTransistionINTO);
                idxOUTOF = find(waveTransistionOUTOF);
                [~,maxIDX] = max(tmpNormalizedSignal);

                % store the peak frame
                TMPpeakFrame(per) = Xi(maxIDX);

                % store the into wave transistion frame
                if isempty(idxINTO)
                    TMPinFrame(per) = 0;
                    TMPvelFrontFrame(per) = 0;
                else
                    % store the wave into
                    TMPinFrame(per) = Xi(idxINTO(1));
                    % store the velocity into
                    TMPvelFrontFrame(per) = tmpVelocitySignal(idxINTO(1));
                end

                % store the into wave transistion frame
                if isempty(idxOUTOF)
                    TMPoutFrame(per) = 0;
                    TMPvelBackFrame(per) = 0;
                else
                    % store the wave outof
                    TMPoutFrame(per) = Xi(idxOUTOF(1));
                    % store the wave velocity outof
                    TMPvelBackFrame(per) = tmpVelocitySignal(idxOUTOF(1));
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % peak frame
            peakFrame(e,:) = TMPpeakFrame;
            % in frame properties
            inFrame(e,:) = TMPinFrame;
            velFrontFrame(e,:) = TMPvelFrontFrame;
            % out frame properties
            outFrame(e,:) = TMPoutFrame;
            velBackFrame(e,:) = TMPvelBackFrame;
        end
        % reshape the inFrame wave transistion
        inFrame = reshape(inFrame,[sz(1:2) numel(waveThresholdPercent)]);
        % reshape the outFrame wave transition
        outFrame = reshape(outFrame,[sz(1:2) numel(waveThresholdPercent)]);
        % reshape the peak frame
        peakFrame = reshape(peakFrame,[sz(1:2) numel(waveThresholdPercent)]);
        % reshape the velocity of the wave front
        velFrontFrame = reshape(velFrontFrame,[sz(1:2) numel(waveThresholdPercent)]);
        % reshape the velocity of the recovery wave
        velBackFrame = reshape(velBackFrame,[sz(1:2) numel(waveThresholdPercent)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        deltaRise = peakFrame - inFrame;
        deltaRecovery = outFrame - peakFrame;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % render t-stacks for inFrame,outFrame,deltaRise,deltaRecovery
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inFrame render
        data = bsxfun(@times,inFrame,plantMask);
        [mM] = getRenderLimits(data);
        renderMap(data,bindVec(plantMask.*STACK(:,:,1)),t_stack_SKIP_frame_vec,[0 1],[1 1],mM,thisOUTPUT);
        % outFrame render
        data = bsxfun(@times,outFrame,plantMask);
        [mM] = getRenderLimits(data);
        renderMap(data,bindVec(plantMask.*STACK(:,:,1)),t_stack_SKIP_frame_vec,[0 1],[1 1],mM,thisOUTPUT);
        % rise time
        data = bsxfun(@times,deltaRise,plantMask);
        [mM] = getRenderLimits(data);
        renderMap(data,bindVec(plantMask.*STACK(:,:,1)),t_stack_SKIP_frame_vec,[0 1],[1 1],mM,thisOUTPUT);
        % recovery time
        data = bsxfun(@times,deltaRecovery,plantMask);
        [mM] = getRenderLimits(data);
        renderMap(data,bindVec(plantMask.*STACK(:,:,1)),t_stack_SKIP_frame_vec,[0 1],[1 1],mM,thisOUTPUT);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out wave analysis to DISK
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for per = 1:numel(waveThresholdPercent)


            tmpFile = [thisOUTPUT '{waveThreshold_' num2str(waveThresholdPercent(per)) '}'];
            csvwrite([tmpFile '{dataType_peakFrame}.csv'],peakFrame(:,:,per));
            csvwrite([tmpFile '{dataType_inFrame}.csv'],inFrame(:,:,per));
            csvwrite([tmpFile '{dataType_outFrame}.csv'],outFrame(:,:,per));
            csvwrite([tmpFile '{dataType_velInWaveFrame}.csv'],velFrontFrame(:,:,per));
            csvwrite([tmpFile '{dataType_velOutWaveFrame}.csv'],velBackFrame(:,:,per));


            propertiesCSV = {};
            ptr = 2;
            propertiesCSV{1,ptr} = 'mean_peak_frame';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_peak_frame';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'mean_in_frame';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_in_frame';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'mean_out_frame';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_out_frame';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'mean_vel_rise';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_vel_rise';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'mean_vel_recover';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_vel_recover';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'mean_rise_time';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_rise_time';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'mean_recover_time';ptr = ptr + 1;
            propertiesCSV{1,ptr} = 'std_recover_time';ptr = ptr + 1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for each mask
            for m = 1:numel(maskFileList)
                % reset the pointer
                ptr = 2;
                % get the path and name for the mask
                [pth,nm] = fileparts(maskFileList{m});
                % index the mask name
                propertiesCSV{m+1,1} = ['{maskFile_' nm];
                % read the r-th mask
                tmp = double(imread(maskFileList{m}));
                tmp = tmp < .5;
                % get the mask idx
                tidx = find(tmp==1);
                % get the stats for peakFrame data
                tmp = peakFrame(:,:,per);
                U_peakFrame = mean(tmp(tidx));
                S_peakFrame = std(tmp(tidx));
                % get the stats for inFrame data
                tmp = inFrame(:,:,per);
                U_inFrame = mean(tmp(tidx));
                S_inFrame = std(tmp(tidx));
                % get the stats for outFrame data
                tmp = outFrame(:,:,per);
                U_outFrame = mean(tmp(tidx));
                S_outFrame = std(tmp(tidx));
                % get the stats for velFrontFrame data
                tmp = velFrontFrame(:,:,per);
                U_velFront = mean(tmp(tidx));
                S_velFront = std(tmp(tidx));
                % get the mean peak velBackFrame data
                tmp = velBackFrame(:,:,per);
                U_velBack = mean(tmp(tidx));
                S_velBack = std(tmp(tidx));
                % get the mean peak deltaRise data
                tmp = deltaRise(:,:,per);
                U_deltaRise = mean(tmp(tidx));
                S_deltaRise = std(tmp(tidx));
                % get the mean peak velBackFrame data
                tmp = deltaRecovery(:,:,per);
                U_deltaRecovery = mean(tmp(tidx));
                S_deltaRecovery = std(tmp(tidx));
                % store in CSV file

                propertiesCSV{m+1,ptr} = U_peakFrame;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_peakFrame;ptr = ptr + 1;

                propertiesCSV{m+1,ptr} = U_inFrame;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_inFrame;ptr = ptr + 1;

                propertiesCSV{m+1,ptr} = U_outFrame;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_outFrame;ptr = ptr + 1;

                propertiesCSV{m+1,ptr} = U_velFront;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_velFront;ptr = ptr + 1;

                propertiesCSV{m+1,ptr} = U_velBack;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_velBack;ptr = ptr + 1;

                propertiesCSV{m+1,ptr} = U_deltaRise;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_deltaRise;ptr = ptr + 1;

                propertiesCSV{m+1,ptr} = U_deltaRecovery;ptr = ptr + 1;
                propertiesCSV{m+1,ptr} = S_deltaRecovery;ptr = ptr + 1;



            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % write out the wave properties in the masks 
            % one sheet for each threshold
            cell2csv([tmpFile '{waveProperties_masks}.csv'],propertiesCSV);
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % track the wave front by isolating the areas
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape SIG
        SIG = reshape(SIG,sz);
        % find the peak time
        sz = size(SIG);
        tmp = reshape(SIG,[prod(sz(1:2)) sz(3)]);
        % loop over all pixels
        for e = 1:size(tmp,1)
            [~,idx] = max(tmp(e,:));
            peakTIME(e) = idx;
        end
        peakTIME = reshape(peakTIME,sz(1:2));
        % masked outside the wave
        m_out_STACK = zeros(size(STACK),'logical');
        % masked inside the wave
        m_in_STACK = zeros(size(STACK),'logical');
        % masked inside sticky wave 
        m_sticky_in_STACK = zeros(size(STACK),'logical');
        %  make the tmp for sticky
        tmp = zeros(size(m_sticky_in_STACK,1),size(m_sticky_in_STACK,2));
        % in number
        inFrame = tmp;
        outFrame = tmp;
        longestIN = tmp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for per = 1:numel(waveThresholdPercent)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make the video writer object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmpFile = [thisOUTPUT '{waveThreshold_' num2str(waveThresholdPercent(per)) '}'];
            videoName = [tmpFile '.avi'];
            v = VideoWriter(videoName);
            open(v)

            for e = 1:globalSKIP:size(STACK,3)
                % mask 'outside' the wave front
                OUTSIDE = SIG(:,:,e) < waveThresholdPercent(per);
                % mask 'inside' the wave front
                INSIDE = SIG(:,:,e) >= waveThresholdPercent(per);
                % masked outside
                maskedOUTSIDE = plantMask.*OUTSIDE;
                % masked inside
                maskedINSIDE = plantMask.*INSIDE;
                % get the stick history
                stickyHISTORY = tmp | INSIDE;
                % current sticky in OR mask inside
                stickyIN = stickyHISTORY | maskedINSIDE;
                % get transitions
                peakMASK = zeros(size(plantMask));
                if e > 1
                    peakMASK = peakTIME.*plantMask == e;
                    [pL1,pL2] = find(peakMASK);

                    %%% record the first time that the wave sweeps over the pixels
                    % transistion into the wave
                    tranINTOwave = m_out_STACK(:,:,e-1) & INSIDE;
                    % index into transistion pixels
                    tidx = find(tranINTOwave);
                    % store the frame at the pixel where the wave starts
                    subIDX = inFrame(tidx) == 0;
                    % store the frame of the transition
                    inFrame(tidx(subIDX)) = e;

                    %%% record the first time that the wave sweeps over the pixels
                    % transistion out of the wave
                    tranOUTOFwave = m_in_STACK(:,:,e-1) & OUTSIDE;
                    % index into transistion pixels
                    tidx = find(tranOUTOFwave);
                    % find those locations where the wave has
                    %   1) never transitioned before
                    %   2) has transistioned into
                    subIDX = outFrame(tidx) == 0 & inFrame(tidx) ~= 0;
                    % store the frame at the pixel where the wave is recovered
                    outFrame(tidx(subIDX)) = e;



                    % make the buffer for the wave front
                    waveFrontBuffer_IN = imdilate(tranINTOwave,strel('disk',1,0));
                    % get the wave front
                    waveFRONT_IN = contourc(double(plantMask.*waveFrontBuffer_IN.*SIG(:,:,e)),[waveThresholdPercent(per) waveThresholdPercent(per)]);

                    % make the buffer for the wave front
                    waveFrontBuffer_OUT = imdilate(tranOUTOFwave,strel('disk',1,0));
                    % get the wave front
                    waveFRONT_OUT = contourc(double(plantMask.*waveFrontBuffer_OUT.*SIG(:,:,e)),[waveThresholdPercent(per) waveThresholdPercent(per)]);


                    %{
                    % show the wave front propogating
                    out = flattenMaskOverlay(double(bindVec(STACK(:,:,e))),tranINTOwave,.75,'r');
                    out = flattenMaskOverlay(out,tranOUTOFwave,.75,'b');
                    imshow(out,[]);
                    title(num2str(e))
                    drawnow
                    %}

                    out = flattenMaskOverlay(double(STACK(:,:,e)),tranINTOwave,.75,'r');
                    out = flattenMaskOverlay(out,tranOUTOFwave,.75,'b');
                    out = flattenMaskOverlay(out,peakMASK,.75,'g');

                    imshow(out,[]);
                    hold on


                    %{
                    % init the pointer
                    cnt = 1;CON = {};
                    % while there are still contours
                    while ~isempty(waveFRONT_IN)
                        % get the contour pointer
                        ptr = waveFRONT_IN(:,1);
                        % get the contour
                        CON{cnt} = waveFRONT_IN(:,2:1+ptr(2));
                        % pop/remove the contour from the queue
                        waveFRONT_IN(:,1:1+ptr(2)) = [];
                        if toRenderVideo
                            % if render is true - render the contour
                            plot(CON{cnt}(1,:),CON{cnt}(2,:),'r');hold on
                        end
                        % increment the contour count
                        cnt = cnt + 1;
                        % store the front for per percent
                        waveFrontContour{e} = CON;
                    end



                    % init the pointer
                    cnt = 1;CON = {};
                    % while there are still contours
                    while ~isempty(waveFRONT_OUT)
                        % get the contour pointer
                        ptr = waveFRONT_OUT(:,1);
                        % get the contour
                        CON{cnt} = waveFRONT_OUT(:,2:1+ptr(2));
                        % pop/remove the contour from the queue
                        waveFRONT_OUT(:,1:1+ptr(2)) = [];
                        if toRenderVideo
                            % if render is true - render the contour
                            plot(CON{cnt}(1,:),CON{cnt}(2,:),'b');hold on
                        end
                        % increment the contour count
                        cnt = cnt + 1;
                        % store the front for per percent
                        waveBackContour{e} = CON;
                    end
                    %}

                    %{
                    if toRenderVideo
                        plot(pL2,pL1,'y.')
                    end
                    %}

                    if toRenderVideo
                        drawnow
                        writeVideo(v,getframe(gca));
                        hold off
                    end

                    %{
                    drawnow
                    hold off
                    %}
                end
                % store the results
                m_out_STACK(:,:,e) = maskedOUTSIDE;
                m_in_STACK(:,:,e) = maskedINSIDE;
                m_sticky_in_STACK(:,:,e) = stickyIN;
                fprintf(['Done finding front in frame ' num2str(e) ':' num2str(size(STACK,3)) '\n']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





            areaD = {};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % init the cell-csv and the rth mask
            for r = 1:numel(maskFileList)
                % make the skip index 
                IDX = (r-1)*3 + 1;
                % get the file parts
                [p,n,ext] = fileparts(maskFileList{r});
                % read the r-th mask
                tmp = double(imread(maskFileList{r}));
                % get the r-th mask
                areaMSK{r} = bindVec(tmp);
                % make the headers for the output csv file
                areaD{1,IDX} = n;
                % make the headers for the 
                areaD{1,IDX+1} = [n '- INSIDE'];
                % make the headers for the 
                areaD{1,IDX+2} = [n '- OUTSIDE'];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % measure each mask for each frame
            % for each mask
            for r = 1:numel(maskFileList)
                % for each frame
                for tm = 1:globalSKIP:size(m_out_STACK,3)
                    %%%%%%%%%%%%%%%%
                    % make the skip index 
                    IDX = (r-1)*3 + 1;
                    %%%%%%%%%%%%%%%%
                    % get filter the inner-wave front by the rth mask
                    tmp = areaMSK{r};
                    % store the area
                    areaD{tm+1,IDX} = sum(tmp(:));
                    %%%%%%%%%%%%%%%%
                    % get filter the inner-wave front by the rth mask
                    tmp = m_in_STACK(:,:,tm).*areaMSK{r};
                    % store the area
                    areaD{tm+1,IDX+1} = sum(tmp(:));
                    %%%%%%%%%%%%%%%%
                    % get filter the inner-wave front by the rth mask
                    tmp = m_out_STACK(:,:,tm).*areaMSK{r};
                    % store the area
                    areaD{tm+1,IDX+2} = sum(tmp(:));
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % write out the cell-csv file
            cell2csv([thisOUTPUT nm '_' baseName '-waveThreshold-' num2str(waveThresholdPercent(per)) '.csv'],areaD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is a "separate" analysis which only draws the contours at the
        % requested levels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the contours and make the video(s) is requested
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the video writer object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        videoName = [oPath nm '_' baseName '.avi'];
        v = VideoWriter(videoName);
        open(v)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:globalSKIP:size(STACK,3)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if render video is true, then show frame
            if toRenderVideo
                imshow(STACK(:,:,e),[]);
                hold on
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % process the contours for each of the percent contour levels
            for p = 1:numel(PER_contour)
                % get the contours at each of the PER_contour levels
                dC = contourc(double(SIG(:,:,e).*plantMask),[PER_contour(p) PER_contour(p)]);
                % init variables for the contours
                cnt = 1;CON = {};
                % while there are still contours
                while ~isempty(dC)
                    % get the contour pointer
                    ptr = dC(:,1);
                    % get the contour
                    CON{cnt} = dC(:,2:1+ptr(2));
                    % pop/remove the contour from the queue
                    dC(:,1:1+ptr(2)) = [];
                    % if render is true - render the contour
                    if toRenderVideo
                        plot(CON{cnt}(1,:),CON{cnt}(2,:),CL{p});hold on
                    end
                    % increment the contour count
                    cnt = cnt + 1;
                end
                % store the contour
                contourStore{p}{e} = CON;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if render is on then save the frame to disk
            if toRenderVideo
                drawnow
                hold off
                writeVideo(v,getframe(gca));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project the wave front(s) as z=stacks
        % over time.  This is helpful to 'see' 
        % the contour evolution rendered on one frame
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project the t-stacks with variable number of skips
        % render to this frame
        toRenderFrame = 1;
        % for each t-SKIP
        for skip = 1:numel(t_stack_SKIP_frame_vec)
            % set the current skip
            SKIP = t_stack_SKIP_frame_vec(skip);
            % for each each percent level
            for p = 1:numel(contourStore)
                % show the first frame
                imshow(STACK(:,:,toRenderFrame),[]);hold on
                % get the number of frames
                NF = numel(1:SKIP:numel(contourStore{p}));
                % get the color map for NF number of frames
                CM = jet(NF);
                % init the contour count
                cnt = 1;
                % for each frame - perform the projection
                for f = 1:SKIP:numel(contourStore{p})
                    % tmp variable for contour
                    tmpCON = contourStore{p}{f};
                    % for each contour - plot in the color from the color map
                    for c = 1:numel(tmpCON)
                        plot(tmpCON{c}(1,:),tmpCON{c}(2,:),'color',CM(cnt,:));
                    end
                    % increment the frame
                    cnt = cnt + 1;
                end
                contourOutPath = [baseTHIS 'CONTOUR_ALGO' filesep];
                mkdir(contourOutPath);
                % save the results of the SKIP t-projections
                saveas(gca,[contourOutPath 'level' num2str(p) '_skip_CONTOUR_ONLY' num2str(SKIP) '_' baseName '.tif']);
            end
            hold off
        end
        close(v)
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








    catch ME
        ME
    end
    
    
    
    
    
end


%{


 %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    find(stickyIN);
    recovery = m_in_STACK;
    sz = size(recovery);
    for e1 = 1:sz(1)
        for e2 = 1:sz(2)
            tmp = recovery(e1,e2,:);
            tmp = bwlarge(tmp);
            r(e1,e2) = sum(tmp);
        end
        e1
    end
    areaD = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % failed attempts at velocity measurements
        SCALES = linspace(1,.05,50);
        filterBank = 5:14;
        VM = [];
        for flt = 1:numel(filterBank)
            %tmpG = imfilter(inFrame,fspecial('gaussian',[61 61],filterBank(flt)),'replicate');
            tmpG = inFrame;
            for scale = 1:numel(SCALES)
                tmp = imresize(tmpG,SCALES(scale));
                tmpMask = imresize(plantMask,SCALES(scale));
                
                [dtdy,dtdx] = gradient(tmp);
                dxdt = dtdx.^-1;
                dydt = dtdy.^-1;
                dxdt(dtdx==0) = 0;
                dydt(dtdy==0) = 0;
                
                velMap = (dydt.^2 + dxdt.^2).^.5;
                
                %{
                [dtdy,dtdx] = gradient(tmp);
                velMap = (dtdx.^-2 + dtdy.^-2).^.5;
                velMap = velMap.*tmpMask;
                %}
                
                imshow(plantMask.*imresize(velMap,size(tmpG)),[]);
                drawnow
                %waitforbuttonpress
                
                VM(flt,scale) = max(velMap(:));
                %plot(log(VM)');
                drawnow
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        % bio-wave calculations
        [dtdy,dtdx] = gradient(inFrame);
        dxdt = dtdx.^-1;
        dydt = dtdy.^-1;
        dxdt(dtdx==0) = 0;
        dydt(dtdy==0) = 0;
        velMap = (dydt.^2 + dxdt.^2).^.5;
        %}
    
        
%}