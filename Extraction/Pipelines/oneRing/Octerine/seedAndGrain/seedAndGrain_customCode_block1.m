function [noReturn] = seedAndGrain_customCode_block1(fileName,output,I,H)
    noReturn = [];
    timingBlock('start','Seed and Grain Image Processing.');
    
    global dataPool;
    
    
    
    I = thawTensor(I);
    H = thawTensor(H);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Get file name parts.');
    if ischar(fileName)
        [p,nm,ext] = fileparts(fileName);
    else
        nm = 'test';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Transform Image to light seed on dark background.');
    if size(I,3) > 1
        H = H(:,:,3);
        % if image is light background with dark seed
        if mean(H(:)) > .5
            H = imcomplement(H);
        end
    else
        if mean(I(:)) > .5
            % if image is light background with dark seed
            H = imcomplement(double(I));
        end
        % make rgb from gray scale
        I = cat(3,I,I,I);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Making mask.');
    % run ostu on value channel - light seed on dark background
    M = H > graythresh(H);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Process mask.');
    % process binary mask
    % remove small objects
    M = bwareaopen(M,50);
    % fill in the object - added for wilson controls of coins
    M = imfill(M,'holes');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Measure with region props.');
    R = regionprops(M,'Area','PixelIdxList','MajorAxisLength','MinorAxisLength','Centroid','Eccentricity');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Running counting method.');
    cidx = count([R.Area]);
    cidx1 = count([R.MajorAxisLength]);
    cidx2 = count([R.MinorAxisLength]);
    fidx = find(cidx==1 & cidx1==1 & cidx2==1);
    singleMask = zeros(size(M));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Create single seed mask.');
    for e = 1:numel(fidx)
        singleMask(R(fidx(e)).PixelIdxList) = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Sampling color channels over all seeds.');
    % for each color channel
    for k = 1:size(I,3)
        % tmp color channel
        tmp = I(:,:,k);
        % for each seed
        for e = 1:numel(fidx)
            cl = double(tmp(R(fidx(e)).PixelIdxList));
            meanRGB(e,k) = mean(cl);
            stdRGB(e,k) = std(cl);
            seedColorSamples(e).colorData(:,k) = cl;
        end
        seedColorSamples(e).Centroid = R(fidx(e)).Centroid;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Make clump mask.');
    % make clump image
    clumpMask = zeros(size(singleMask));
    clidx = find(cidx~=1);
    totalCount = cidx;
    % remove clumps greater than 100 seeds
    totalCount(totalCount > 100) = 0;
    % get clump count
    totalCount = sum(totalCount);
    % make clump mask
    for e = 1:numel(clidx)
        clumpMask(R(clidx(e)).PixelIdxList) = 1;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Make overlay single(red) and clump (green) overlay.');
    out = flattenMaskOverlay(I,logical(singleMask),.6,'r');
    out = flattenMaskOverlay(out,logical(clumpMask),.6,'g');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isdeployed()
        numSeedsToSample = numel(seedColorSamples);
    else
        numSeedsToSample = min(10,numel(seedColorSamples));
    end
    timingBlock('note',['Writing out whole color sheets per seed.....' num2str(numSeedsToSample) '.']);
    centroidData = [];
    allColors = [];
    for seed = 1:numSeedsToSample
        tmp = seedColorSamples(seed).colorData;
        tmp = tmp(randperm(size(tmp,1)),:);
        samp = min(2000,size(tmp,1));
        allColors = cat(1,allColors,tmp(1:samp,:));


        %seedFileName_color = ['{seedNumber_' num2str(seed) '}{feature_color}.csv'];
        clear fileStruct
        fileStruct.fileName = nm;
        fileStruct.seedNumber = num2str(seed);
        fileStruct.feature = 'color';
        fileStruct.ext = 'csv';
        seedFileName_color = outPort.formatFileName(fileStruct);

        spoolData('user','local',seedFileName_color,seedColorSamples(seed).colorData,{});
        centroidData = [centroidData;seedColorSamples(seed).Centroid];
    end

    %seedFileName_centroid = ['{seedNumber_all}{feature_centroid}.csv'];
    clear fileStruct
    fileStruct.fileName = nm;
    fileStruct.seedNumber = 'all';
    fileStruct.feature = 'centroid';
    fileStruct.ext = 'csv';
    seedFileName_centroid = outPort.formatFileName(fileStruct);
    spoolData('user','local',seedFileName_centroid,centroidData,{});
      

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note',['Write phenotype data.']);
    phenotypeData = [[R(fidx).Area]' [R(fidx).MajorAxisLength]' [R(fidx).MinorAxisLength]' [R(fidx).Eccentricity]' meanRGB stdRGB];
    

    %seedPhenotypeData_fileName = [nm '_phenotypeData.csv'];
    clear fileStruct
    fileStruct.fileName = nm;
    fileStruct.ext = 'csv';
    fileStruct.type = 'phenotypeData';
    seedPhenotypeData_fileName = outPort.formatFileName(fileStruct);
    spoolData('user','local',seedPhenotypeData_fileName,phenotypeData,{});

     
    
    
    U_data = mean(phenotypeData,1);
    STD_data = std(phenotypeData,1,1);

    
    
    
    [U,E,L] = PCA_FIT_FULLws(allColors,1,false);
    alpha = 1;
    LS = alpha*linspace(-L.^.5,L.^.5,5);
    colorSquare = [];
    for s = 1:numel(LS)
        tmpC = [];
        for k = 1:3
            tmpC(:,:,k) = (U(k)+E(k)*LS(s))*ones(101,101);
            eyeColor(s,k) = U(k)+E(k)*LS(s);
        end
        colorSquare = [colorSquare;tmpC];
    end
    colorSquare = colorSquare/255;
    
    
    %colorPatch_fileName = [filesep 'colorPatch.jpg'];
    clear fileStruct
    fileStruct.ext = 'jpg';
    fileStruct.fileName = nm;
    fileStruct.type = 'colorPatch';
    colorPatch_fileName = outPort.formatFileName(fileStruct);
    spoolData('user','local',colorPatch_fileName,colorSquare,{});
    
    %colorPatchData_fileName = [filesep 'colorPatch.csv'];
    clear fileStruct
    fileStruct.ext = 'csv';
    fileStruct.fileName = nm;
    fileStruct.type = 'colorPatch';
    colorPatchData_fileName = outPort.formatFileName(fileStruct);
    spoolData('user','local',colorPatchData_fileName,eyeColor,{});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('note','Write image data.');
    %returnImage_fileName = [nm '_return.jpg'];
    clear fileStruct
    fileStruct.ext = 'jpg';
    fileStruct.fileName = nm;
    fileStruct.type = 'return';
    returnImage_fileName = outPort.formatFileName(fileStruct);
    spoolData('user','local',returnImage_fileName,out,{});
    
    
    
    
    timingBlock('start','Creating JSON file');
    tmpDoc = [];
    timingBlock('note','Adding shape information to JSON file.');
    tmpDoc = generatePhenotypeNode(tmpDoc,nm,'FileName','FileName');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(1),'AverageArea','AverageArea');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(1),'StdArea','StdArea');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(2),'MajorAxisLength','MajorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(2),'StdMajorAxisLength','StdMajorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(3),'MinorAxisLength','MinorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(3),'StdMinorAxisLength','StdMinorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(4),'AverageEccentricity','AverageEccentricity');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(4),'StdEccentricity','StdEccentricity');
    
    timingBlock('note','Adding color information to JSON file.');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(5),'AverageMeanRed','AverageMeanRed');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(5),'StdMeanRed','StdMeanRed');
    
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(6),'AverageMeanGreen','AverageMeanGreen');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(6),'StdMeanGreen','StdMeanGreen');
    
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(7),'AverageMeanBlue','AverageMeanBlue');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(7),'StdMeanBlue','StdMeanBlue');
    
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(8),'AverageStdRed','AverageStdRed');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(8),'StdStdRed','StdStdRed');
    
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(9),'AverageStdGreen','AverageStdGreen');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(9),'StdStdGreen','StdStdGreen');
    
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(10),'AverageStdBlue','AverageStdBlue');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(10),'StdStdBlue','StdStdBlue');
    
    for s = 1:size(eyeColor,1)
        tmpDoc  = generatePhenotypeNode(tmpDoc,eyeColor(s,1),...
            ['Red_colorRange_step' num2str(s)],['Red_colorRange_step' num2str(s)]);
        tmpDoc  = generatePhenotypeNode(tmpDoc,eyeColor(s,2),...
            ['Green_colorRange_step' num2str(s)],['Green_colorRange_step' num2str(s)]);
        tmpDoc  = generatePhenotypeNode(tmpDoc,eyeColor(s,3),...
            ['Blue_colorRange_step' num2str(s)],['Blue_colorRange_step' num2str(s)]);
    end
    
    timingBlock('note','Adding total count information to JSON file.');
    tmpDoc = generatePhenotypeNode(tmpDoc,totalCount,'totalCount','totalCount');

    


    JSON_string = savejson('seedDoc',tmpDoc);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % save JSON string
    %%%%%%%%%%%%%%%%%%%%%%%%%
    jsonOUTfile = [nm '_jdoc.json'];
    fileID = fopen(jsonOUTfile,'w');
    fprintf(fileID,strrep(JSON_string,'\/','\\/'));
    fclose(fileID);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    timingBlock('stop');
    
    
    
    timingBlock('stop');
end