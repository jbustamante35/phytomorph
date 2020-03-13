function [] = extractModelDataFromWaveMovie(fileName,oPath)

    %[pth,nm,ext] = fileparts(fileName);
    %oPath = [pth filesep 'extraction' filesep];
    %oPath = [basePath filesep nm filesep];
    mkdir(oPath);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter for spatial gradient
    filt_SIZE = [11 11];
    filt_STD = 3;
    filt_spatial = fspecial('gaussian',filt_SIZE,filt_STD);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filter for load
    filt_load = fspecial('gaussian',[41 41],11);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % down sample number
    DN = 10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % wave threshold
    waveThreshold = .3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of points along rise and recover
    TMI = 150;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reporting the timing of each step
    tmg = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load and stack data
    % can be either a CZI of TIF file
   
    STACK = czi_tif_loader(fileName,filt_load);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STACK = double(STACK);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalize the raw stack
    [SIG] = normalizeRawStack(STACK);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make whole plant mask
    SS = std(STACK,1,3);
    % make the mean signal of the image
    U = mean(STACK,3);
    % make the plant mask based on the standard dev
    plantMask = bindVec(SS) > graythresh(bindVec(SS))*.2;
    plantMask = bwlarge(plantMask);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create the gradient data and filter
    [d1,d2] = gradient(SIG);
    for tm = 1:size(SIG,3)
        d1(:,:,tm) = imfilter(d1(:,:,tm),filt_spatial,'replicate');
        d2(:,:,tm) = imfilter(d2(:,:,tm),filt_spatial,'replicate');
    end
    dN = (d1.^2 + d2.^2).^-.5;
    d1 = d1.*dN;
    d2 = d2.*dN;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sz = size(SIG);
    SIG = reshape(SIG,[prod(sz(1:2)) sz(3)]);
    d1 = reshape(d1,[prod(sz(1:2)) sz(3)]);
    d2 = reshape(d2,[prod(sz(1:2)) sz(3)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oSIG = reshape(SIG,sz);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [s1,s2] = ndgrid(1:size(plantMask,1),1:size(plantMask,2));
    ds = mod(s1,DN) == 0 & mod(s2,DN) == 0;
    didx = find(ds);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plantMask = bwlarge(plantMask);
    pidx = find(plantMask);
    pidx = intersect(pidx,didx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp = false;
    [n2,n1] = ndgrid(linspace(-20,20,50),linspace(-20,20,50));
    domainSZ = size(n1);
    domain = [n1(:) n2(:) ones(size(n1(:))) ones(size(n2(:)))]';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Q1 = zeros([domainSZ TMI numel(pidx)]);
    %Q2 = zeros([domainSZ TMI numel(pidx)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subSIG = SIG(pidx,:);
    subd1 = d1(pidx,:);
    subd2 = d2(pidx,:);
    [subP(:,1),subP(:,2)] = ind2sub(sz(1:2),pidx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
    kp = [];
    parfor e = 1:numel(pidx)
        try
            TTtm = clock;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm  = clock;
            sigStack = [subSIG(e,:);subd1(e,:);subd2(e,:)];
            etm = etime(clock,tm);
            if tmg;fprintf(['stack data signals :' num2str(etm) '\n']);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm  = clock;
            %[p(1),p(2)] = ind2sub(sz(1:2),pidx(e));
            p = subP(e,:);
            etm = etime(clock,tm);
            if tmg;fprintf(['translate p-idx to p-xy :' num2str(etm) '\n']);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm  = clock;
            [riseDomain,recoverDomain] = generateBinaryDomain(sigStack,waveThreshold);
            etm = etime(clock,tm);
            if tmg;fprintf(['generate binary domains  :' num2str(etm) '\n']);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm  = clock;
            [rise_affineSequence] = generateAffineSequence(sigStack,p,riseDomain);
            [recover_affineSequence] = generateAffineSequence(sigStack,p,recoverDomain);
            etm = etime(clock,tm);
            if tmg;fprintf(['generate affine sequences :' num2str(etm) '\n']);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm = clock;
            [rise_affineSequence] = rescaleData(rise_affineSequence,TMI);
            [recover_affineSequence] = rescaleData(recover_affineSequence,TMI);
            etm = etime(clock,tm);
            if tmg;fprintf(['rescale affine data :' num2str(etm) '\n']);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm  = clock;
            riseFrames = oSIG(:,:,(riseDomain));
            recoverFrames = oSIG(:,:,(recoverDomain));
            etm = etime(clock,tm);
            fprintf(['index frames2 :' num2str(etm) '\n'])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm = clock;
            [rise_data_out] = coreSample(oSIG,domain,rise_affineSequence,domainSZ,disp);
            [recover_data_out] = coreSample(oSIG,domain,recover_affineSequence,domainSZ,disp);
            etm = etime(clock,tm);
            if tmg;fprintf(['sample core data :' num2str(etm) '\n']);end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tm = clock;
            [rise_data_out] = rescaleData(rise_data_out,200);
            [recover_data_out] = rescaleData(recover_data_out,200);
            etm = etime(clock,tm);
            fprintf(['rescale core data :' num2str(etm) '\n'])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
            writeToDisk(oPath,rise_data_out,recover_data_out,rise_affineSequence,recover_affineSequence,p)
            
            %Q1(:,:,:,e) = rise_data_out;
            %Q2(:,:,:,e) = recover_data_out;
            etm = etime(clock,TTtm);
            fprintf(['total time :' num2str(etm) '\n'])
            kp(e) = true;
        catch
            kp(e) = false;
        end
    end
    
    
    
    matFile = [oPath '{image_plantMask}.mat'];
    save(matFile,'plantMask','kp');
    
end
%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dig for the raw data for extracting the wave front data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FilePath = '/mnt/snapper/nate/Analysis/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the base path for the extracted data
basePath = '/mnt/snapper/nate/AminoAcidWave/Extraction/';
mkdir(basePath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run through the data - create fail list, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
failList = {};
preDoneList = {};
cnt2 = 1;
cnt1 = 1;
for e = 1:numel(FileList)
    try
        % get the file parts
        [pth,nm,ext] = fileparts(FileList{e});
        fidx = strfind(pth,filesep);
        % create the oPath for the extracted data
        oPath = [basePath pth((fidx(end)+1):end) filesep];
        % if the path exists then do not run
        if ~exist(oPath)
            extractModelDataFromWaveMovie(FileList{e},oPath);
        else
            preDoneList{cnt1} = FileList{e};
            cnt1 = cnt1 + 1;
        end
    catch ME
        % catch the fail list
        failList{cnt2} = FileList{e};
        cnt2 = cnt2 + 1;
    end
end


%}