%% load first set
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/fieldEarData/ear_images/';
FileList = {};
FileExt = {'jpg'};
FileList = fdig(FilePath,FileList,FileExt,1);
%% load second set
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/fieldEarData/fromOverSeas/';
FileList = {};
FileExt = {'jpg'};
FileList = fdig(FilePath,FileList,FileExt,1);
%% loop view images
for e = 1:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%% algorithm
close all
areaThresh = 4000;
oPath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/fieldEarData/fromOverSeas/return1/';
mmkdir(oPath);

for e = 1:numel(FileList)
    fileName = FileList{e};
    processFieldEarImage(fileName,oPath,areaThresh);
end
