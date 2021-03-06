%% list out the contents for the return folder
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/results/{projectName_TIPS}/{functionName_main}/{createDate_2020-04-28-22-02-33}/';
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/results/{projectName_TIPS}/{functionName_main_ver2}/{createDate_2020-04-29-16-04-48}/';
FileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG','png'};
FileList = dig(FilePath,FileList,FileExt,1);
%% look through the pile of data
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/tassel/UIUC_2019_Tassel_images/';
FileList = {};
FileExt = {'jpg'};
rawFileList = dig(FilePath,FileList,FileExt,1);
%% find the data which is done with compute
close all
oPath = '/mnt/spaldingdata/nate/octerineDataStore/results/{projectName_TIPS}/{functionName_main}/{createDate_2020-04-28-22-02-33}/';
maskTemplateName = [oPath '{originalName_<FILE>}{maskNumber_<MASK>}.tif'];
processedTemplateName = [oPath '{originalName_<FILE>}{maskNumber_<MASK>}{fileType_processedImage}.png'];
cnt = 1;
disp = false;
clear data notLoad
cntN = 1;
for e = 1:numel(rawFileList)
    
    bool = [];
    
    fileName = rawFileList{e};
    [p,n,ex] = fileparts(fileName);
    
    maskName = {};
    parfor m = 1:4
        maskName{m} = strrep(strrep(maskTemplateName,'<FILE>',n),'<MASK>',num2str(m));
        bool(m) = exist(maskName{m},'file');
    end
    processedName = {};
    parfor m = 1:4
        processedName{m} = strrep(strrep(processedTemplateName,'<FILE>',n),'<MASK>',num2str(m));
        %bool = [bool exist(processedName{m},'file')];
        bool(4+m) = exist(processedName{m},'file');
    end
    
    
    if all(bool)

        area = [];
        MN = [];
        for m = 1:4
            M = imread(maskName{m});
            area(1,m) = sum(M(:));
            M = logical(M);
            M = imclearborder(M);
            area(2,m) = sum(M(:));
        end

        data(cnt).fileName = fileName;
        data(cnt).area = area;
        data(cnt).processedName = processedName;
        data(cnt).maskName = maskName;
        cnt = cnt + 1;
        
    else
        notLoad(cntN).fileName = rawFileList;
        notLoad(cntN).bool = bool';
        cntN = cntN + 1;
    end
    
    
    
    if (any(bool) & ~all(bool))
        break
    end
    
    if disp
        if all(bool)
            I = imread(rawFileList{e});

            MN = [];
            for m = 1:4
                M = imread(maskName{m});
                out = flattenMaskOverlay(I,M);
                MN = cat(4,MN,out);
            end

            PN = [];
            for m = 1:4
                M = imread(processedName{m});
                M = permute(M,[2 1 3]);
                M = imresize(M,[size(I,1) size(I,2)]);
                PN = cat(4,PN,M);
            end

            MN = cat(4,MN,PN);

            
            
            MN = montage(MN,'Size', [2 4]);
            drawnow
        end
    end
    fprintf(['Done with:' num2str(e) ':' num2str(numel(rawFileList)) '\n'])
end
%%
for e = 1:numel(notLoad)
    notLoad(e).bool
end
%% scan for data from data structure
Area = [];
for e = 1:numel(data)
    Area(e) = sum(abs(diff(data(e).area,1,1)));
end
fidx = find(Area ~= 0);
clear failList
for e = 1:numel(fidx)
    failList(e) = data(fidx(e));
end
%%
fileName = '/mnt/spaldingdata/nate/octerineDataStore/tassel/UIUC_2019_Tassel_images/DDPSC_02010-1-90_0256.jpg';
func = project.getFunction('main');
func.operate(fileName,makeTempLocation())