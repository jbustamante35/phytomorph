function [] = measureArabidopsisSeed(fileName,oPath)
    fileList = {};
    % make directory
    fprintf(['Make output directory.\n']);
    mkdir(oPath);
    fprintf(['Read Image.\n']);
    if ischar(fileName)
        I = imread(fileName);
    else
        I = fileName;
    end
    % remove the "alpha-channel" - added Oct 25. 2019
    if size(I,3) > 3
        I(:,:,4:end) = [];
    end
    
    
    fprintf(['Transform RGB-->HSV.\n']);
    if size(I,3) > 1
         H = rgb2hsv(I);
         H = H(:,:,3);
         if mean(H(:)) > .5
            H = imcomplement(H);
         end
    else
         H = imcomplement(double(I))/255+1;
         I = cat(3,I,I,I);
    end
    fprintf(['Threshold image.\n']);
    fprintf(['Image size:' num2str(size(H)) '\n']);
    % run ostu on value channel
    M = H > graythresh(H);
    % remove small objects
    M = bwareaopen(M,50);
    % fill in the object - added for wilson controls of coins
    M = imfill(M,'holes');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Measure with region props.\n']);
    R = regionprops(M,'Area','PixelIdxList','MajorAxisLength','MinorAxisLength','Centroid','Eccentricity');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Running counting method.\n']);
    cidx = count([R.Area]);
    cidx1 = count([R.MajorAxisLength]);
    cidx2 = count([R.MinorAxisLength]);
    fidx = find(cidx==1 & cidx1==1 & cidx2==1);
    singleMask = zeros(size(M));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Running counting method.\n']);
    for e = 1:numel(fidx)
        singleMask(R(fidx(e)).PixelIdxList) = 1;
    end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Sampling color channels over all seeds.\n']);
    meanRGB = zeros(numel(fidx),size(I,3));
    stdRGB = meanRGB;
    for k = 1:size(I,3)
        tmp = I(:,:,k);
        for e = 1:numel(fidx)
            cl = double(tmp(R(fidx(e)).PixelIdxList));
            meanRGB(e,k) = mean(cl);
            stdRGB(e,k) = std(cl);
            seedColorSamples(e).colorData(:,k) = cl;
            seedColorSamples(e).Centroid = R(fidx(e)).Centroid;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make clump image
    clumpMask = zeros(size(singleMask));
    clidx = find(cidx~=1);
    totalCount = cidx;
    totalCount(totalCount > 100) = 0;
    totalCount = sum(totalCount);
    for e = 1:numel(clidx)
        clumpMask(R(clidx(e)).PixelIdxList) = 1;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Make overlay.\n']);
    out = flattenMaskOverlay(I,logical(singleMask),.6,'r');
    out = flattenMaskOverlay(out,logical(clumpMask),.6,'g');


    %{
    image(out);
    hold on
    for e = 1:numel(clidx)
        text(R(clidx(e)).Centroid(1)-100,R(clidx(e)).Centroid(2)-100,num2str(cidx(clidx(e))),'Background','w');
        plot([R(clidx(e)).Centroid(1)-100 R(clidx(e)).Centroid(1)],[R(clidx(e)).Centroid(2)-100 R(clidx(e)).Centroid(2)],'y')
    end
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Get file parts.\n']);
    if ischar(fileName)
        [p,nm,ext] = fileparts(fileName);
    else
        nm = 'test';
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Write image data.\n']);
    fileList{end+1} = [oPath nm '_return.jpg'];
    imwrite(out,fileList{end});

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Added all seed color for Craine Brothers
    % Due to algorithm taking a long time to run
    % this feature is being turned off for now - Dec, 12, 2019
    if isdeployed()
        numSeedsToSample = numel(seedColorSamples);
    else
        numSeedsToSample = min(10,numel(seedColorSamples));
    end
    fprintf(['Writing out whole color sheets per seed.....' num2str(numSeedsToSample) '\n'])
    seedFileName_centroid = [oPath '{seedNumber_all}{feature_centroid}.csv'];
    centroidData = [];
    allColors = [];
    for seed = 1:numSeedsToSample
        tmp = seedColorSamples(seed).colorData;
        tmp = tmp(randperm(size(tmp,1)),:);
        samp = min(2000,size(tmp,1));

        allColors = cat(1,allColors,tmp(1:samp,:));
        
        
        %seedFileName_color = [oPath '{seedNumber_' num2str(seed) '}{feature_color}.csv'];
        %csvwrite(seedFileName_color,seedColorSamples(seed).colorData);
        
        centroidData = [centroidData;seedColorSamples(seed).Centroid];
    end
    csvwrite(seedFileName_centroid,centroidData);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Write phenotype data.\n']);
    phenotypeData = [[R(fidx).Area]' [R(fidx).MajorAxisLength]' [R(fidx).MinorAxisLength]' [R(fidx).Eccentricity]' meanRGB stdRGB];
    fileList{end+1} = [oPath nm '_phenotypeData.csv'];
    csvwrite(fileList{end},phenotypeData);

    
    U_data = mean(phenotypeData,1);
    STD_data = std(phenotypeData,1,1);

    
    
    
    
    
    [U,E,L] = PCA_FIT_FULLws(allColors,1);
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
    imwrite(colorSquare,[oPath filesep 'colorPatch.jpg']);
    csvwrite([oPath filesep 'colorPatch.csv'],eyeColor);
    
    
    
    
    
    tmpDoc = [];
    tmpDoc = generatePhenotypeNode(tmpDoc,nm,'FileName','FileName');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(1),'AverageArea','AverageArea');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(1),'StdArea','StdArea');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(2),'MajorAxisLength','MajorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(2),'StdMajorAxisLength','StdMajorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(3),'MinorAxisLength','MinorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(3),'StdMinorAxisLength','StdMinorAxisLength');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(4),'AverageEccentricity','AverageEccentricity');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(4),'StdEccentricity','StdEccentricity');
    
    
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
    
    
    tmpDoc = generatePhenotypeNode(tmpDoc,totalCount,'totalCount','totalCount');

    


    JSON_string = savejson('seedDoc',tmpDoc);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % save JSON string
    %%%%%%%%%%%%%%%%%%%%%%%%%
    fileList{end+1} = [oPath nm '_jdoc.json'];
    fileID = fopen(fileList{end},'w');
    fprintf(fileID,strrep(JSON_string,'\/','\\/'));
    fclose(fileID);
    %%%%%%%%%%%%%%%%%%%%%%%%%

end

%{
    fileName = '/home/nate/10_26_trials/GS_test_COL_4800dpi_rgb028.tif';
    oPath = '/home/nate/10_26_trials/';
    measureArabidopsisSeed(fileName,oPath);

    I = imread('/home/nate/10_26_trials/GS_test_COL_4800dpi_rgb028_return.jpg');

    FilePath = '/home/nate/10_26_trials/Calibration Images/';
    oPath = '/home/nate/10_26_trials/Calibration Images_return/';
    mkdir(oPath);
    FileList = {};
    FileExt = {'tiff','TIF'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    for e = 1:numel(FileList)
        measureArabidopsisSeed(FileList{e},oPath);
    end

    fileName = '/iplant/home/turnersarahd/Clement_drought_scans/Normal/035255243056';
    measureArabidopsisSeed(fileName,'');

    fileName = '~/1010.tif'';
    measureArabidopsisSeed(fileName,'');

    fileName = '/home/nate/Blue_wheat.jpg';
    measureArabidopsisSeed(fileName,'');



    fileName = '/mnt/snapper/nate/Col_bottomplate_3sheetswhitepaper.jpg';
    fileName = '/mnt/snapper/nate/Col_bottomplate_3papers.tif';
    fileName = '/mnt/snapper/nate/Bob.tif';
    fileName = '/mnt/snapper/nate/btx623.tiff';
    measureArabidopsisSeed(fileName,'./output/');
    
%}


