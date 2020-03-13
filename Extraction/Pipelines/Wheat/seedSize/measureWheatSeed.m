function [] = measureWheatSeed(fileName,oPath)
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
    fprintf(['Transform RGB-->HSV.\n']);
    if size(I,3) > 1
         H = rgb2hsv(I);
         H = H(:,:,3);
    else
         H = imcomplement(double(I))/255+1;
         I = cat(3,I,I,I);
    end
    fprintf(['Threshold image.\n']);
    fprintf(['Image size:' num2str(size(H)) '\n']);
    
    Lab = rgb2lab(I);
    isColor = mean(abs(Lab(:,:,2:3)),3);
    
    % run ostu on value channel
    M = H > graythresh(H);
    % remove small objects
    M = bwareaopen(M,50);
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
    for k = 1:size(I,3)
        tmp = I(:,:,k);
        for e = 1:numel(fidx)
            cl = double(tmp(R(fidx(e)).PixelIdxList));
            meanRGB(e,k) = mean(cl);
            stdRGB(e,k) = std(cl);
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
    fprintf(['Write phenotype data.\n']);
    phenotypeData = [[R(fidx).Area]' [R(fidx).MajorAxisLength]' [R(fidx).MinorAxisLength]' [R(fidx).Eccentricity]' meanRGB stdRGB];
    fileList{end+1} = [oPath nm '_phenotypeData.csv'];
    csvwrite(fileList{end},phenotypeData);

    
    U_data = mean(phenotypeData,1);
    STD_data = std(phenotypeData,1,1);

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
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(6),'AverageStdGreen','AverageStdGreen');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(7),'AverageMeanBlue','AverageMeanBlue');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(7),'AverageStdBlue','AverageStdBlue');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(8),'AverageStdRed','AverageStdRed');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(8),'StdStdRed','StdStdRed');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(9),'AverageStdGreen','AverageStdGreen');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(9),'StdStdGreen','StdStdGreen');
    tmpDoc = generatePhenotypeNode(tmpDoc,U_data(10),'AverageStdBlue','AverageStdBlue');
    tmpDoc = generatePhenotypeNode(tmpDoc,STD_data(10),'StdStdBlue','StdStdBlue');
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
    measureWheatSeed(fileName,'');


    
%}


