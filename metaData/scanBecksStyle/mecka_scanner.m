%% dig json files for the cobs
cjFilePath = '/mnt/tetra/nate/becksData/WCIC/processedResults/cobData/';
cjFileList = {};
cjFileExt = {'json'};
cjFileList = gdig(cjFilePath,cjFileList,cjFileExt,1);
%% dig json files for the ears
ejFilePath = '/mnt/tetra/nate/becksData/WCIC/processedResults/earData/';
ejFileList = {};
ejFileExt = {'json'};
ejFileList = gdig(ejFilePath,ejFileList,ejFileExt,1);
%% dig json files for the kernels
kjFilePath = '/mnt/tetra/nate/becksData/processedResults/kernelData/';
kjFilePath = '/mnt/tetra/nate/becksData/processedResults/kernelData/BecksKernelData-2019-03-06-20-33-29.4/';
kjFilePath = '/mnt/tetra/nate/becksData/WCIC/processedResults/kernelData/BecksKernelData-2019-03-06-20-33-29.4/';
% old one
%kjFilePath = '/mnt/tetra/nate/becksData/WCIC/processedResults/kernelData/phytoMorph_Image_Phenomics_Tool_Kit_analysis1-2019-02-06-12-56-25.4/';
kjFileList = {};
kjFileExt = {'json'};
kjFileList = gdig(kjFilePath,kjFileList,kjFileExt,1);
%% find the fail numbers for cobs
expectedNumber = 394;
analysisNumberCob = [];
for e = 1:numel(cjFileList)
    fidx = strfind(cjFileList{e},'analysis-');
    fidx2 = strfind(cjFileList{e},filesep);
    fidx2(fidx2 < fidx(1)) = [];
    analysisNumberCob(e) = str2num(cjFileList{e}((fidx(1) + 9):(fidx2(1)-1)));
end
errorCobs = setdiff(1:expectedNumber,analysisNumberCob)
numel(errorCobs)


%% find the fail numbers for ears
expectedNumber = 394*2;
analysisNumberEar = [];
for e = 1:numel(ejFileList)
    fidx = strfind(ejFileList{e},'analysis-');
    fidx2 = strfind(ejFileList{e},filesep);
    fidx2(fidx2 < fidx(1)) = [];
    analysisNumberEar(e) = str2num(ejFileList{e}((fidx(1) + 9):(fidx2(1)-1)));
end
errorEars = setdiff(1:expectedNumber,analysisNumberEar);
%% find the fail numbers for kernels
expectedNumber = 394;
analysisNumberCob = [];
for e = 1:numel(kjFileList)
    fidx = strfind(kjFileList{e},'analysis-');
    fidx2 = strfind(kjFileList{e},filesep);
    fidx2(fidx2 < fidx(1)) = [];
    analysisNumberKernel(e) = str2num(kjFileList{e}((fidx(1) + 9):(fidx2(1)-1)));
end
errorKernels = setdiff(1:expectedNumber,analysisNumberKernel);
numel(errorKernels)
errorKernels
%% get cob list
cFilePath = '/mnt/tetra/nate/becksData/cobData/';
cFileList = {};
cFileExt = {'tif'};
cFileList = gdig(cFilePath,cFileList,cFileExt,1);
errFileName = '/mnt/tetra/nate/becksData/processedResults/cobData/BecksCobData-2019-02-05-13-28-11.6/analysis-XXX/logs/condor-stdout-0';
errorSearchName = {};
other = [];
for e = 1:numel(errorCobs)
    try
        eName = strrep(errFileName,'XXX',num2str(errorCobs(e)));
        errLog = readtext(eName);
        for l = 1:size(errLog,1)
            if ~isempty(errLog{l,1})
                if ischar(errLog{l,1})
                    if contains(errLog{l,1},'End ticket strip:')
                        fidx = strcmp(errLog{l,1},':');
                        tmpName = errLog{l,1}((18):end);
                        if ~isempty(tmpName)
                            errorSearchName{end+1} = tmpName;
                        end
                    end
                end
            end
        end
    catch
        other = [other ;e];
    end
end
errorSearchName = unique(errorSearchName);
%%
close all
for s = 1:numel(errorSearchName)
    
    toSearch = '14178600_imageData_cob.tif';
    toSearch = '14182011_imageData_cob.tif';
    toSearch = errorSearchName{s};
    for e = 1:numel(cFileList)
        if (contains(cFileList{e},toSearch))
            fidx(e) = true;
        else
            fidx(e) = false;
        end
    end
    I = imread(cFileList{find(fidx)});
    imshow(I,[]);
    drawnow
    
end
%% get ear list
eFilePath = '/mnt/tetra/nate/becksData/earData/';
eFileList = {};
eFileExt = {'tif'};
eFileList = gdig(eFilePath,eFileList,eFileExt,1);
toSearch = '14178600_imageData_cob.tif';
toSearch = '14182011_imageData_cob.tif';
for e = 1:numel(eFileList)
    if (contains(eFileList{e},toSearch))
        fidx(e) = true;
    else
        fidx(e) = false;
    end
end
I = imread(eFileList{find(fidx)});
imshow(I,[]);
%% get kernel list
kFilePath = '/mnt/tetra/nate/becksData/kernelData/';
kFileList = {};
kFileExt = {'tif'};
kFileList = gdig(kFilePath,kFileList,kFileExt,1);
toSearch = '14182008_imageData_kernel.tif';
toSearch = '15725897_imageData_kernel.tif';
for e = 1:numel(kFileList)
    if (contains(kFileList{e},toSearch))
        fidx(e) = true;
    else
        fidx(e) = false;
    end
end
I = imread(kFileList{find(fidx)});
imshow(I,[]);
%% auto search for kernel fails
kFilePath = '/mnt/tetra/nate/becksData/cobData/';
kFileList = {};
kFileExt = {'tif'};
kFileList = gdig(kFilePath,kFileList,kFileExt,1);
errFileName = '/mnt/tetra/nate/becksData/processedResults/kernelData/BecksKernelData-2019-03-06-20-33-29.4/analysis-XXX/logs/condor-stdout-0';
errorSearchName = {};
other = [];
for e = 1:numel(errorKernels)
    try
        eName = strrep(errFileName,'XXX',num2str(errorKernels(e)));
        errLog = readtext(eName);
        for l = 1:size(errLog,1)
            if ~isempty(errLog{l,1})
                if ischar(errLog{l,1})
                    if contains(errLog{l,1},'FileName:')
                        fidx = strcmp(errLog{l,1},':');
                        tmpName = errLog{l,1}((10):end);
                        if ~isempty(tmpName)
                            errorSearchName{end+1} = tmpName;
                        end
                    end
                end
            end
        end
    catch
        other = [other ;e];
    end
end
errorSearchName = unique(errorSearchName);
%%
close all
for s = 1:numel(errorSearchName)
    
    toSearch = '14178600_imageData_cob.tif';
    toSearch = '14182011_imageData_cob.tif';
    toSearch = errorSearchName{s};
    for e = 1:numel(kFileList)
        if (contains(kFileList{e},toSearch))
            fidx(e) = true;
        else
            fidx(e) = false;
        end
    end
    I = imread(cFileList{find(fidx)});
    imshow(I,[]);
    drawnow
    
end

%%
runFileList = kFileList(find(fidx));
    mecka('k',runFileList{1},3,'/home/nate/Downloads/',[],1,1,1200,1);
%%


